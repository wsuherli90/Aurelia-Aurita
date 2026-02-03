package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

class AxionFieldUnit extends Component {
  val io = new Bundle {
     val ricci_scalar = slave Stream(RealType())
     val det_g        = slave Stream(RealType())
     
     val axion_out    = master Stream(RealType())
  }
  
  val join = StreamJoin.arg(io.ricci_scalar, io.det_g)
  
  val s1 = join.translateWith {
      val r = join.payload._1
      val d = join.payload._2
      
      val alpha = r * d 
      alpha
  }.s2mPipe()
  
  io.axion_out << s1
}
