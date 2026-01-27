package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

/**
 * Active Stress Unit
 * Computes T_ij = mu * g_ij
 * Input: Stream[Bundle]
 * Output: Stream[Matrix3]
 */
class ActiveStressUnit extends Component {
  val io = new Bundle {
    val input_valid = in Bool()
    val output_ready = in Bool()
    
    val mu_in     = in(RealType())
    val g_metric  = in(Matrix3())
    
    val output = master Stream(Matrix3())
  }
  
  case class StressInput() extends Bundle {
     val mu = RealType()
     val g  = Matrix3()
  }
  
  val input_stream = Stream(StressInput())
  input_stream.valid := io.input_valid
  input_stream.payload.mu := io.mu_in
  input_stream.payload.g  := io.g_metric

  val s1 = input_stream.translateWith {
      val i = input_stream.payload
      val t = Matrix3()
      for(r <- 0 until 3; c <- 0 until 3) {
          t(r)(c) := i.mu * i.g(r)(c)
      }
      t
  }.s2mPipe()
  
  io.output << s1
}
