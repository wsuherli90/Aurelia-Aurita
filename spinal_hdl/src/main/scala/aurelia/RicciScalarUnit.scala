package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

/**
 * Ricci Scalar Unit
 * Computes R = g^ij * R_ij
 * Input: Stream[Bundle]
 * Output: Stream[SFix]
 */
class RicciScalarUnit extends Component {
  val io = new Bundle {
    val input_valid = in Bool()
    val output_ready = in Bool()
    
    // Inputs (Bundled logic externally, but here we define explicit ports)
    val ricci_tensor = in(Matrix3()) // R_ij
    val g_inv        = in(Matrix3()) // g^ij
    
    val output = master Stream(RealType())
  }
  
  case class ScalarInput() extends Bundle {
     val ricci = Matrix3()
     val g_inv = Matrix3()
  }
  
  val input_stream = Stream(ScalarInput())
  input_stream.valid := io.input_valid
  input_stream.payload.ricci := io.ricci_tensor
  input_stream.payload.g_inv := io.g_inv
  
  // Pipeline Stage 1: Parallel Multiplication
  val s1 = input_stream.translateWith {
      val i = input_stream.payload
      val products = Vec(RealType(), 9)
      for(r <- 0 until 3; c <- 0 until 3) {
          products(r*3 + c) := i.g_inv(r)(c) * i.ricci(r)(c)
      }
      products
  }.s2mPipe()
  
  // Pipeline Stage 2: Summation (Adder Tree)
  val s2 = s1.translateWith {
      val products = s1.payload
      
      // Reduce logic
      // We can use a registered tree if needed, but s2mPipe covers the path delay
      val sum = products.reduce(_ + _)
      sum
  }.s2mPipe()
  
  io.output << s2
}

