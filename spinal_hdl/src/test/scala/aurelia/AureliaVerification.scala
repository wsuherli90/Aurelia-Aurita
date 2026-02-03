package aurelia

import spinal.core._
import spinal.core.sim._
import spinal.lib.bus.amba4.axis._

object AureliaVerification {
  def main(args: Array[String]): Unit = {
    SimConfig.withWave.compile(new AureliaAccelerator()).doSim { dut =>
      dut.clockDomain.forkStimulus(period = 10)
      
      dut.io.s_axis_data.valid #= false
      dut.io.s_axis_data.last #= false
      dut.io.m_axis_data.ready #= true 
      
      dut.io.grid_dim(0) #= 2
      dut.io.grid_dim(1) #= 2
      dut.io.grid_dim(2) #= 2
      
      val TOTAL_VOXELS = 8
      
      dut.clockDomain.waitSampling(10)
      
      println(s"Sending $TOTAL_VOXELS Beats of Data (Full Frame)...")
      
      dut.clockDomain.waitSampling()
      dut.io.s_axis_data.valid #= true
      
      for(i <- 0 until TOTAL_VOXELS) {
         dut.io.s_axis_data.payload.data #= BigInt(i) * BigInt("10000000000000") 
         dut.io.s_axis_data.last #= (i == TOTAL_VOXELS - 1)
         
         dut.clockDomain.waitSampling()
         while(!dut.io.s_axis_data.ready.toBoolean) {
            dut.clockDomain.waitSampling()
         }
      }
      
      dut.io.s_axis_data.valid #= false
      dut.io.s_axis_data.last #= false
      
      println("Data Sent. Waiting for Output...")
      
      var received_count = 0
      var received_last = false
      
      for(_ <- 0 until 100) {
         if (dut.io.m_axis_data.valid.toBoolean) {
            received_count += 1
            println(s"Received Beat $received_count: Data=${dut.io.m_axis_data.payload.data.toBigInt.toString(16)} Last=${dut.io.m_axis_data.payload.last.toBoolean}")
            
            if (dut.io.m_axis_data.payload.last.toBoolean) {
                received_last = true
                assert(received_count == TOTAL_VOXELS, s"TLAST asserted early/late! Count=$received_count Expected=$TOTAL_VOXELS")
            }
         }
         dut.clockDomain.waitSampling()
      }
      
      assert(received_count == TOTAL_VOXELS, s"Did not receive all data! Got $received_count")
      assert(received_last, "TLAST was never asserted!")
      println("Verification SUCCESS: AXI4 Stream & TLAST Verified.")
    }
  }
}
