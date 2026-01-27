package aurelia

import spinal.core._
import spinal.sim._
import spinal.core.sim._
import aurelia.AureliaTypes._

object SimRicciFlow {
  def main(args: Array[String]): Unit = {
    // Compile the Top Level Accelerator
    SimConfig.withWave.compile(new AureliaAccelerator()).doSim { dut =>
      dut.clockDomain.forkStimulus(period = 10)

      // Initialize Inputs
      dut.io.host_data_in.valid #= false
      dut.io.host_data_out.ready #= true // Always ready (Sink)
      dut.io.y_vec_in.valid #= false
      
      // Configure Grid Dimensions (Small for Sim)
      // NX=10, NY=10, NZ=5 -> 500 Voxels
      dut.io.grid_config(0) #= 10
      dut.io.grid_config(1) #= 10
      dut.io.grid_config(2) #= 5
      
      val total_voxels = 10 * 10 * 5 

      dut.clockDomain.waitSampling(10)

      var sent_count = 0
      var received_count = 0
      
      // Monitor Output
      fork {
         while(received_count < total_voxels) {
             dut.io.host_data_out.ready #= true
             dut.clockDomain.waitSampling()
             
             if(dut.io.host_data_out.valid.toBoolean) {
                 received_count += 1
                 if(received_count % 50 == 0) {
                     println(s"[SIM] Received Voxel $received_count / $total_voxels")
                 }
                 // Simple sanity check: Output shouldn't be all zero if logic works?
                 // Or it might be, if stress is zero.
                 // We just check valid flow for now.
             }
         }
         println(s"[SIM] SUCCESS! Processed $received_count voxels.")
         simSuccess()
      }

      // Drive Input (Host Stream + Y Vector)
      // SFix Range: 16 exp, -16 exp. 
      // 1.0 = 2^16 = 65536
      val ONE_FIXED = 65536L
      
      while(sent_count < total_voxels) {
          dut.io.host_data_in.valid #= true
          dut.io.y_vec_in.valid #= true
          
          // Fill Symmetric Matrix with Identity-like pattern
          // 00=1, 11=1, 22=1. Others=0.
          // Indices: 0,3,5 are diagonals.
          for(i <- 0 until 6) {
             if(i == 0 || i == 3 || i == 5) {
                dut.io.host_data_in.payload(i).raw #= ONE_FIXED 
             } else {
                dut.io.host_data_in.payload(i).raw #= 0
             }
          }
          
          // Fill Y-Vector (Physics placeholder)
          dut.io.y_vec_in.payload(0).raw #= 0
          dut.io.y_vec_in.payload(1).raw #= 0
          dut.io.y_vec_in.payload(2).raw #= 0
          
          dut.clockDomain.waitSampling()
          
          // Handshake logic
          if(dut.io.host_data_in.ready.toBoolean) {
             sent_count += 1
          }
      }
      
      dut.io.host_data_in.valid #= false
      dut.io.y_vec_in.valid #= false
      
      // Timeout watchdog
      dut.clockDomain.waitSampling(3000)
      println("[SIM] TIMEOUT! Pipeline stalled or latency too high.")
      simFailure()
    }
  }
}
