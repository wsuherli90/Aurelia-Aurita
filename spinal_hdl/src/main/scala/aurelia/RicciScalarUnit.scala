package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

class RicciScalarUnit extends Component {
  val io = new Bundle {
    val input_valid = in Bool()
    val output_ready = in Bool()
    
    val ricci_tensor = in(Matrix3()) 
    val g_inv        = in(Matrix3()) 
    
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
  
  val s1 = input_stream.translateWith {
      val i = input_stream.payload
      val products = Vec(RealType(), 9)
      for(r <- 0 until 3; c <- 0 until 3) {
          products(r*3 + c) := i.g_inv(r)(c) * i.ricci(r)(c)
      }
      products
  }.s2mPipe()
  
  val s2 = s1.translateWith {
      val products = s1.payload
      
      val sum = products.reduce(_ + _)
      sum
  }.s2mPipe()
  
  io.output << s2
}
