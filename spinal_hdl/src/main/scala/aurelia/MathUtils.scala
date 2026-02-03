package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

object MathUtils {
  class ReciprocalUnit extends Component {
     val io = new Bundle {
       val input  = slave Stream(RealType())
       val output = master Stream(RealType())
     }
     
     // Algorithm: Pipelined Newton-Raphson Division with Manual Stall Handling
     // Note: "Manual Stall" logic (RegNextWhen) is used here to ensure precise control 
     // over validity propagation in this complex pipeline.
     // Future Refactor: Wrap in Stream.pipelined for better SpinalHDL style.
     
     val stall = !io.output.ready
     
     val sign = io.input.payload.raw.msb
     val abs_input = io.input.payload.abs
     val lzc = CountLeadingZeros(abs_input.raw.asBits)
     
     // S1
     val s1_valid = RegNextWhen(io.input.valid, !stall) init(False)
     val s1_val   = RegNextWhen(io.input.payload, !stall)
     val s1_lzc   = RegNextWhen(lzc, !stall)
     val s1_sign  = RegNextWhen(sign, !stall)
     
     // S2 (Pre-Scale & LUT/Linear approx)
     // To avoid complex LUT, use x0 = 2.91 - 2x approx for [0.5, 1]
     // Norm val:
     val s2_valid = RegNextWhen(s1_valid, !stall) init(False)
     val norm_raw = s1_val.abs.raw |<< s1_lzc
     val norm_x   = RealType(); norm_x.raw := norm_raw
     // Linear approx x0 = 2.823 - 1.882 * x (Chebyshev for 1/x on [0.5, 1])
     val s2_x0    = RealConst(2.8235) - (norm_x * RealConst(1.8823))
     val s2_lzc   = RegNextWhen(s1_lzc, !stall)
     val s2_sign  = RegNextWhen(s1_sign, !stall)
     
     // S3: NR 1
     // x1 = x0 * (2 - x * x0)
     val s3_val = RegNextWhen(norm_x, !stall) // Need x for NR
     val s3_term = RegNextWhen(norm_x * s2_x0, !stall)
     val s3_x0   = RegNextWhen(s2_x0, !stall)
     
     val s3_valid = RegNextWhen(s2_valid, !stall) init(False)
     val s3_params = (RegNextWhen(s2_lzc, !stall), RegNextWhen(s2_sign, !stall))

     // S4: NR 1 result
     // sub = 2 - term
     // x1 = x0 * sub
     val s4_sub = RealConst(2.0) - s3_term
     val s4_x1  = s3_x0 * s4_sub
     val s4_valid = RegNextWhen(s3_valid, !stall) init(False)
     val s4_params = (RegNextWhen(s3_params._1, !stall), RegNextWhen(s3_params._2, !stall))
     val s4_x_orig = RegNextWhen(s3_val, !stall)

     // S5: NR 2
     val s5_valid = RegNextWhen(s4_valid, !stall) init(False)
     val s5_term  = RegNextWhen(s4_x_orig * s4_x1, !stall)
     val s5_x1    = RegNextWhen(s4_x1, !stall)
     val s5_params = (RegNextWhen(s4_params._1, !stall), RegNextWhen(s4_params._2, !stall))
     
     // S6: NR 2 Result
     val s6_sub  = RealConst(2.0) - s5_term
     val s6_x2   = s5_x1 * s6_sub
     val s6_valid = RegNextWhen(s5_valid, !stall) init(False)
     val s6_params = (RegNextWhen(s5_params._1, !stall), RegNextWhen(s5_params._2, !stall))
     
     // Output Check and Re-scaling
     // 1/(x * 2^-k) = 1/x * 2^k.
     // input was x_orig. norm = x_orig * 2^lzc.
     
     val final_val = s6_x2 |<< s6_params._1
     val final_signed = s6_params._2 ? -final_val | final_val
     
     io.output.valid   := s6_valid
     io.output.payload := final_signed
     io.input.ready    := !stall 
  }
  
  // Wrapper for implicit Flow-like usage where backpressure is known to be handled
  // or the consumer is guaranteed to accept data (e.g. valid-only interface).
  def reciprocal(input: SFix, valid: Bool): (SFix, Bool) = {
     val u = new ReciprocalUnit()
     u.io.input.payload := input
     u.io.input.valid   := valid
     u.io.output.ready  := True 
     (u.io.output.payload, u.io.output.valid)
  }

      // Input Setup: x * 2^16
      // Input width = 32 (16.16). We shift left by 16 to treat it as Q32.32 before sqrt?
      // No. To get Q16.16 output, we need Q32.32 input equivalent scale.
      // Sqrt( I * 2^-16 ) = Sqrt(I) * 2^-8. We want 2^-16 out.
      // So we need Sqrt(X) * 2^-16.
      // Sqrt( I * 2^-32 ) = Sqrt(I) * 2^-16.
      // So input 'I' must be scaled by 2^-32 relative to integer 1.
      // Input is Q16.16 => I * 2^-16. We need I * 2^-16 * 2^16 = I (raw).
      // If we shift input left by 16 bits, we get scale 2^-32.
      // Example: 1.0 = 0x00010000. Shift << 16 = 0x100000000.
      // Sqrt(0x100000000) = 0x10000 = 1.0 in Q16.16. Correct.
      // Input width 32 + 16 shift = 48 bits. Sqrt(48 bits) -> 24 bits.
      
      val PIPELINE_STAGES = 24
      
      // Initial state
      case class SqrtState(rem: UInt, root: UInt, valid: Bool)
      
      val stages = Vec(SqrtState(0,0,False), PIPELINE_STAGES + 1)
      
      // Stage 0
      val x_extended = io.input.payload.raw.asUInt.resize(48) |<< 16
      stages(0).rem   := x_extended
      stages(0).root  := 0
      stages(0).valid := io.input.valid
      io.input.ready  := True // Simplified pipeline: assuming always ready or handled at end?
      // WARNING: Stream must handle backpressure. If we don't support backpressure, we lose compliant AXI.
      // But implementing backpressure on a fine-grained arithmetic pipeline is expensive.
      // Standard solution: Backpressure stops the whole pipeline.
      // We will use a `halt` signal.
      
      val halt_pipeline = !io.output.ready
      
      for(i <- 0 until PIPELINE_STAGES) {
         val current = stages(i)
         val next    = stages(i+1)
         
         // Logic for bit i (from MSB 23 down to 0)
         val bit     = 23 - i 
         val shift   = bit * 2
         
         // We look at the top bits of Remainder corresponding to current window?
         // Standard Digit Recurrence:
         // R_{i+1} = 2*R_i - d_i * (2Q_i + d_i * 2^-i?)
         // Let's use the standard integer algorithm adapted.
         // Current Root Q. We check if (2Q + 1) <= (R >> shift).
         
         // Actually, R is shifted left by 2 each step in standard algo.
         // Keep R static and shift masks? No, shift R is better.
         
         // Let's restart logic cleanly for "Restoring Square Root"
         // At each step, we bring down 2 bits of the radicant.
         // But we have the full radicant available at start.
      }
  }
  
  // Impl: Simple unrolled restoring square root
  // Note: One stage per iteration unrolled. For very high frequency (>300MHz),
  // this logic path might need further breaking (retiming) depending on FPGA fabric speed.
  class FixedSqrtPipelined extends Component {
      val io = new Bundle {
          val input  = slave Stream(SFix(16 exp, -16 exp))
          val output = master Stream(SFix(16 exp, -16 exp))
      }
      
      // We process 2 bits per stage. Total 48 input bits -> 24 stages?
      // 48 bits input. Output 24 bits.
      // Let's assume input is positive.
      
      val din = io.input.payload.raw.asUInt.resize(48) << 16
      
      // Pipeline registers
      case class Stage(rem: UInt, root: UInt)
      val pipes = Vec(Reg(Stage(0, 0)), 25) 
      val valids = Vec(RegInit(False), 25)
      
      // Pipeline stall control
      val stall = !io.output.ready
      
      when(!stall) {
         valids(0) := io.input.valid
         pipes(0).rem  := din
         pipes(0).root := 0
         
         for(i <- 0 until 24) {
            val rem  = pipes(i).rem
            val root = pipes(i).root
            
            // Look at top 2 bits of Remainder?
            // Actually, we usually shift Rem at each step.
            // Let's shift Rem Left by 2 at each step.
            // And pull in 0s? No, logic is:
            // Rem matches the input bits.
            // Let's use the subtraction based method:
            // Bit k (from 23 down to 0).
            // Test Val = (root << 1 | 1) << k?
            // This is confusing. 
            // Standard approach: 
            // val test = (root << 2) + 1. If (rem >> (46 - 2*i)) >= test?
            
            // Correct approach:
            // At start, rem = input.
            // loop i from 0 to 23:
            //   shift = 46 - 2*i
            //   group = (rem >> shift) // but we are consuming bits!
            
            // BETTER: Shift-based updates.
            // R_{next} = (R_{curr} << 2) | next_2_bits_of_input
            // But doing that across pipeline requires carrying input.
            // EASIER: Remainder holds current residue.
            // At step 0: Rem = Input >> 46?
            
            // Let's rely on simple large logic per stage if needed or just correct math.
            // Given the task, I will instantiate a simple iterative logic but unrolled.
            
            // Stage i: determining bit (23-i) of Root.
            // Let `rem` be the working remainder.
            // Let `root` be the accumulated root so far.
            
            // In iteration 0, we deal with bits 47..46.
            // BUT: We need to pass the FULL original input down? YES or shift it out.
            // Let's pass `rem` as the stream of input bits remaining.
            
            // Actually, `pipes` definition above is insufficient.
            // Redefine:
            pipes(i+1).rem  := pipes(i).rem // Placeholder
            valids(i+1) := valids(i)
         }
      }
      
      io.input.ready := !stall
      
      // Real Implementation of logic inside the loop:
      // We will use a fully unrolled function on the signal to generate the next state.
      
      def step(stage: Int, r: UInt, q: UInt): (UInt, UInt) = {
          // r: current remainder (contains all bits, effectively shifted)
          // q: current root
          
          // We are solving for bit `b = 23 - stage`.
          // Value to subtract: T = (q << 2) | 1.
          // Where is T aligned? T is aligned with the top of R?
          // Standard algo: 
          // R = (R << 2) | next_bits ... if we ingest.
          // If R is full width:
          // Check if R >= (T << (2*b)).
          // If yes, R -= (T << (2*b)); q |= (1 << b).
          
          val b = 23 - stage
          val t = (q << 1) | 1
          val t_shifted = t << b // This shift is incorrect scaling.
          
          // Let's use "Restoring Division" style square root.
          // bit b (from 23 to 0).
          // We test subtraction of (2*Q + 2^b)*2^b?
          // (Q + 2^b)^2 = Q^2 + 2*Q*2^b + 2^2b.
          // Subtracted Q^2 already.
          // Need to subtract 2*Q*2^b + 2^2b = (2*Q + 2^b) * 2^b.
          // Q is already shifted?
          
          // Simplified:
          // Q is accumulated root above bit b.
          // Test term: ( (Q << 1) | 1 ) << b . 
          // Note: Q here is `root` so far, treated as integer `root >> (b+1)`?
          // No, Q is the raw integer building up.
          
          // Correct Loop:
          // for(b = 23; b >= 0; b--)
          //   term = (root << 1) | (1 << b) ? No.
          //   val test = (root | (1<<b))
          //   if (test*test <= input) -> too expensive (multipliers).
          
          // Back to Subtraction method:
          // rem is maintained.
          // b = 23. test = (1 << 23).
          // term = (root << 1) | (1 << b). 
          // term = term << b.
          // if rem >= term: rem -= term, root |= (1<<b).
          
          // Let's apply this.
          // `root` is accumulating. `rem` is reducing.
          val shift = b
          val term_base = (q << 1) | (U(1) << shift) // wait, q has bits > b set? yes.
          // No, q is `root` so far.
          // q has bits 23..(b+1) set.
          // We want to add bit 'b'.
          // (q + 2^b)^2 - q^2 = 2*q*2^b + (2^b)^2 = (2*q + 2^b) * 2^b
          // term = ( (q << 1) + (1<<b) ) << b
          
          val term = ((q << 1) | (U(1) << b)) << b
          
          val diff = r - term
          // Check valid subtraction (r >= term)
          val can_sub = r >= term
          
          val next_r = can_sub ? diff | r
          val next_q = can_sub ? (q | (U(1) << b)) | q
          
          (next_r, next_q)
      }
      
      // Apply unrolled logic to pipeline registers
      when(!stall) {
          for(i <- 0 until 24) {
             val (nr, nq) = step(i, pipes(i).rem, pipes(i).root)
             pipes(i+1).rem  := nr
             pipes(i+1).root := nq
          }
      }
      
      io.output.valid := valids(24)
      io.output.payload.raw := pipes(24).root.asSInt.resize(32)
  }
}
}
