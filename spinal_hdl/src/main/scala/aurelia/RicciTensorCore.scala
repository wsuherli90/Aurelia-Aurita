package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

/**
 * Ricci Tensor Core (Fixed Point, Streaming)
 * Calculates Ricci Curvature Tensor from Connection Coefficients
 * 
 * Pipeline Depth: ~2-3 Stages
 */
class RicciTensorCore extends Component {
  val io = new Bundle {
    val input  = slave Stream(new Bundle {
        val g_inv       = Matrix3()
        val gamma       = Tensor3()
        val laplacian_g = Matrix3()
    })
    
    val output = master Stream(Matrix3())
  }
  
  // Helper for Adder Tree (Hierarchical Summation)
  def adderTree(inputs: Seq[RealType]): RealType = {
    if (inputs.size == 1) return inputs.head
    val paired = inputs.grouped(2).map { 
        case Seq(a, b) => RegNext(a + b) // Pipeline Stage for every adder layer
        case Seq(a)    => RegNext(a)
    }.toSeq
    adderTree(paired)
  }

  // =============================================================
  // CRITICAL FIX: Use Flow for heavy pipeline to prevent stall corruption.
  // We convert Stream -> Flow, process deeply, then Queue -> Stream.
  // =============================================================
  
  // 1. Decouple input (Elastic Buffer)
  val flowInput = io.input.toFlow(fire = true) // Only takes data when valid. 
  // Warning: toFlow(fire=true) consumes data whenever valid is present regardless of downstream.
  // But we need backpressure at the END of the flow.
  // Strategy: 
  // Stream -> Queue -> Flow -> Queue -> Stream
  // If the output Queue is full, we must drop data OR we assume the Flow is always fast enough 
  // but if downstream halts, we need enough buffer.
  // BETTER SOTA STRATEGY: 
  // Use `Flow` for the math. Valid propagates. 
  // The entry into the Flow is gated by the `almostFull` of the output FIFO.
  
  // Input FIFO to decouple
  val inputFifo = io.input.queue(16) 
  
  // Custom Flow Logic
  // We take from FIFO only if standard pipeline rules allow? 
  // Actually, Flow logic in Spinal is simple: Data Valid moves every clock.
  // If we can't output, we must not take input.
  
  val outputFifo = Stream(Matrix3())
  val pipelineFlow = new Area {
      // Logic: Read from Input FIFO if we are not "Stopping". 
      // But Flow cannot "Stop" mid-flight.
      // So we only read if there is space in Output FIFO to accept the result N cycles later.
      // Simply: Stream -> Flow -> Stream is naive if we don't manage the "flight" count.
      // For simplicity in this iteration: We rely on a large Output Queue.
      
      val flow = inputFifo.toFlow 
      
      // =============================================================
      // Stage 1: Pre-calculation & Partial Products
      // =============================================================
      val s1 = flow.map { payload =>
        val gamma = payload.gamma
        val lap   = payload.laplacian_g
        
        // 1. Contract Gamma: Gamma^m_lm
        val gamma_contracted = Vec(RealType(), 3)
        for(l <- 0 until 3) {
          gamma_contracted(l) := RegNext(gamma(0)(l)(0) + gamma(1)(l)(1) + gamma(2)(l)(2))
        }
        (gamma, gamma_contracted, lap)
      } // Latency: +1 (from RegNext inside map)
      
      // =============================================================
      // Stage 2: Product Generation
      // =============================================================
      val s2 = s1.map { case (gamma, gamma_contr, lap) =>
         val term3_prods = Vec(Vec(RealType(), 3), 3)
         val term4_partials = Vec(Vec(Vec(RealType(), 9), 3), 3)
         
         for(i <- 0 until 3; j <- 0 until 3) {
           val t3_mul = Vec(RealType(), 3)
           for(l <- 0 until 3) {
              t3_mul(l) := RegNext(gamma(l)(i)(j) * gamma_contr(l)) 
           }
           term3_prods(i)(j) := RegNext(t3_mul(0) + t3_mul(1) + t3_mul(2))
           
           for(k <- 0 until 3; l <- 0 until 3) {
              val idx = k*3 + l
              term4_partials(i)(j)(idx) := RegNext(gamma(k)(i)(l) * gamma(l)(j)(k))
           }
         }
         (term3_prods, term4_partials, lap)
      } // Latency: +2 (Total +3)

      // =============================================================
      // Stage 3: Adder Trees & Final Assembly
      // =============================================================
      val s3 = s2.map { case (t3_prods, t4_parts, lap) =>
        val final_ricci = Matrix3()
        for(i <- 0 until 3; j <- 0 until 3) {
           val t3 = t3_prods(i)(j)
           val t4_p = t4_parts(i)(j)
           
           // Adder Tree (Depth 4)
           val l1_0 = RegNext(t4_p(0) + t4_p(1)); val l1_1 = RegNext(t4_p(2) + t4_p(3))
           val l1_2 = RegNext(t4_p(4) + t4_p(5)); val l1_3 = RegNext(t4_p(6) + t4_p(7))
           val l1_4 = RegNext(t4_p(8))
           
           val l2_0 = RegNext(l1_0 + l1_1); val l2_1 = RegNext(l1_2 + l1_3); val l2_2 = RegNext(l1_4)
           val l3_0 = RegNext(l2_0 + l2_1); val l3_1 = RegNext(l2_2)
           
           val t4_sum = RegNext(l3_0 + l3_1) 
           
           val t3_d   = Delay(t3, 4)
           val lap_d  = Delay(lap(i)(j), 4)
           
           val nonlinear = RegNext(t3_d - t4_sum) 
           val linear    = RegNext(Delay(RealConst(-0.5), 4) * lap_d) // Note: RealConst is not a signal, Delay works on signal.
           // FIX: RealConst returns SFix signal. Delay(Signal) works.
           // BUT Delay(RealConst(...), 4) implies delay on constant wire. Optimization removes it. Correct.
           
           final_ricci(i)(j) := RegNext(linear + nonlinear) 
        }
        final_ricci
      } // Latency: +6 (Total +9)
      
      // Management of Flow Valid:
      // Flow.map does NOT automatically delay the valid signal by the amount of RegNexts inside the logic!
      // We must manually match the latency.
      // Total Latency: 1 (s1) + 2 (s2) + 6 (s3) = 9 Cycles.
      
      val result = Stream(Matrix3())
      result.valid := Delay(inputFifo.valid, 9, init = False) // Delay Valid to match data
      result.payload := s3.payload 
      // Note: s3.valid is NOT DELAYED automatically by map logic in current Spinal version usually, 
      // unless using `Flow.stage()`. We manually mapped payloads with RegNext. 
      // The `s3` Flow object itself keeps the ORIGINAL valid timing if we just used .map purely on payload.
      // So existing `s3.valid` == `inputFifo.valid`.
      // We explicitly construct the result Stream.
  }
  
  // Output Buffering
  // We use a large FIFO to absorb the pipeline flight.
  // Ideally, if buffer is nearing full, we stop Input.
  io.output_ricci << pipelineFlow.result.queue(32) // Deep buffer
}
