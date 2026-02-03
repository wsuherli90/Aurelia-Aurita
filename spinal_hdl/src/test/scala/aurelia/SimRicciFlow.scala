package aurelia

import spinal.core._
import spinal.sim._
import spinal.core.sim._
import aurelia.AureliaTypes._

object SimRicciFlow {
  import scala.util.Random
  
  // Golden Model Function (Double Precision)
  def goldenRicciFlow(g_center: Array[Double], ricci_val: Array[Double], dt: Double): Array[Double] = {
      val res = new Array[Double](6)
      // g_new = g_old - 2 * dt * ricci
      // Simplified check: We assume we verify the integration step given a Ricci tensor.
      // But we need to verify the Ricci tensor calculation itself?
      // That requires full neighborhood. 
      // For this test, we verify the Identity Preservation Property:
      // If Input is Flat Space (Identity), Ricci should be 0. Output should be Identity.
      
      // We also verify AXI robustness by random throttling.
      for(i <- 0 until 6) {
         res(i) = g_center(i) // Expect no change for flat space
      }
      res
  }

  def main(args: Array[String]): Unit = {
    SimConfig.withWave.compile(new AureliaAccelerator()).doSim { dut =>
      dut.clockDomain.forkStimulus(period = 10)

      dut.io.s_axis_data.valid #= false
      dut.io.m_axis_data.ready #= true 
      dut.io.grid_dim(0) #= 10
      dut.io.grid_dim(1) #= 10
      dut.io.grid_dim(2) #= 5
      
      val total_voxels = 10 * 10 * 5 

      // Robustness: Fork a process that randomly toggles m_axis_data.ready (Backpressure)
      fork {
         while(true) {
            dut.io.m_axis_data.ready #= Random.nextBoolean()
            dut.clockDomain.waitSampling()
         }
      }

      var sent_count = 0
      var received_count = 0
      var error_count = 0
      
      // Monitor Process
      fork {
         val ONE_DOUBLE = 1.0
         val TOLERANCE  = 0.001
         
         while(received_count < total_voxels) {
             dut.clockDomain.waitSampling()
             
             if(dut.io.m_axis_data.valid.toBoolean && dut.io.m_axis_data.ready.toBoolean) {
                 val payload = dut.io.m_axis_data.payload.data.toBigInt
                 // Unpack 6x 32-bit values
                 // We sent Matrix Identity. Expected: Identity.
                 
                 // SFix 16.16: 1.0 = 65536 (0x10000)
                 // Check diagonal elements
                 val diag_0 = (payload >> 0) & 0xFFFFFFFFL
                 val diag_1 = (payload >> 32) & 0xFFFFFFFFL
                 val diag_2 = (payload >> 64) & 0xFFFFFFFFL
                 
                 val expected = 65536
                 
                 // Simple assertion for flat space logic
                 if( (diag_0 - expected).abs > 10 || (diag_1 - expected).abs > 10 || (diag_2 - expected).abs > 10 ) {
                     println(s"[FAIL] Voxel $received_count: Expected Identity, Got $diag_0, $diag_1, $diag_2")
                     error_count += 1
                 }
                 
                 received_count += 1
                 if(received_count % 50 == 0) println(s"[SIM] Progress: $received_count/$total_voxels")
             }
         }
         
         if(error_count == 0) {
             println(s"[SIM] SUCCESS! Verified Flat Space preservation under AXI backpressure.")
         } else {
             println(s"[SIM] FAILURE! Found $error_count errors.")
             simFailure()
         }
         simSuccess()
      }

      val ONE_FIXED = 65536L // 1.0 in 16.16
      
      // Helper to pack 512 bits
      def packInput(g: Array[Long], y: Array[Long]): BigInt = {
         var res = BigInt(0)
         for(i <- 0 until 6) res |= (BigInt(g(i)) << (i*32))
         for(i <- 0 until 3) res |= (BigInt(y(i)) << ((6+i)*32))
         res
      }
      
      while(sent_count < total_voxels) {
          dut.io.s_axis_data.valid #= true
          
          // Send Identity Matrix + Zero Vector
          val g_in = Array(ONE_FIXED, ONE_FIXED, ONE_FIXED, 0L, 0L, 0L)
          val y_in = Array(0L, 0L, 0L)
          
          dut.io.s_axis_data.payload.data #= packInput(g_in, y_in)
          
          dut.clockDomain.waitSampling()
          
          if(dut.io.s_axis_data.ready.toBoolean) {
             sent_count += 1
          }
      }
      
      dut.io.s_axis_data.valid #= false
      
      dut.clockDomain.waitSampling(5000)
      println("[SIM] TIMEOUT! Pipeline stalled.")
      simFailure()
    }
  }
}
