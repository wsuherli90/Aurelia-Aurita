package aurelia

import spinal.core._
import spinal.sim._
import spinal.core.sim._
import AureliaTypes._

import org.scalatest.funsuite.AnyFunSuite

class ZeroSkippingTest extends AnyFunSuite {
  test("ZeroSkippingScheduler should skip flat spacetime") {
    SimConfig.withWave.compile(new ZeroSkippingScheduler()).doSim { dut =>
      dut.clockDomain.forkStimulus(period = 10)
      
      dut.io.threshold.raw #= 10 
      val threshold_raw = (0.1 * (1<<16)).toInt
      dut.io.threshold.raw #= threshold_raw
      
      dut.io.input_window.valid #= false
      dut.clockDomain.waitSampling(5)
      
      val ident_val = (1.0 * (1<<16)).toInt
      
      def setFlat() = {
          for(z<-0 until 3; y<-0 until 3; x<-0 until 3) {
             for(r<-0 until 3; c<-0 until 3) {
                if(r==c) dut.io.input_window.payload(z)(y)(x)(r)(c).raw #= ident_val
                else     dut.io.input_window.payload(z)(y)(x)(r)(c).raw #= 0
             }
          }
      }
      
      dut.io.input_window.valid #= true
      setFlat()
      
      dut.clockDomain.waitSampling(1) 
      
      var skipped = false
      for(i <- 0 until 10) {
         dut.clockDomain.waitSampling()
         if(dut.io.skip_cmd.valid.toBoolean && dut.io.skip_cmd.payload.toBoolean) {
             skipped = true
         }
      }
      assert(skipped, "Failed to skip flat spacetime")
      
      val perturbed_val = ident_val + (0.2 * (1<<16)).toInt
      dut.io.input_window.payload(1)(1)(1)(0)(0).raw #= perturbed_val
      
      dut.clockDomain.waitSampling(1)
      
      var working = false
      for(i <- 0 until 10) {
         dut.clockDomain.waitSampling()
         if(dut.io.skip_cmd.valid.toBoolean && !dut.io.skip_cmd.payload.toBoolean) {
             working = true
         }
      }
       assert(working, "Failed to process curved spacetime (False Positive Skip)")
    }
  }
}
