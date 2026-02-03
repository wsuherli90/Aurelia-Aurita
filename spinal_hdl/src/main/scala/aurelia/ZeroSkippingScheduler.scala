package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

class ZeroSkippingScheduler extends Component {
  val io = new Bundle {
    val input_window = slave Stream(Window3D())
    val skip_cmd     = master Stream(Bool()) 
    
    val threshold    = in(RealType()) 
  }

  val s1 = io.input_window.translateWith {
    val w = io.input_window.payload
    val center = w(1)(1)(1)
    
    val neighbors = Vec(
      w(0)(1)(1), w(2)(1)(1), 
      w(1)(0)(1), w(1)(2)(1), 
      w(1)(1)(0), w(1)(1)(2)  
    )
    
    val diffs = neighbors.map { n =>
       val sum_sq = RealType()
       
       val partials = Vec(RealType(), 9)
       for(r <- 0 until 3; c <- 0 until 3) {
         val d = center(r)(c) - n(r)(c)
         partials(r*3+c) := d * d
       }
       
       val norm = partials.reduceBalancedTree(_ + _)
       norm
    }
    
    val total_variance = diffs.reduceBalancedTree(_ + _)
    total_variance
  }.s2mPipe() 

  val s2 = s1.translateWith {
    val variance = s1.payload
    val is_flat = variance < io.threshold 
    is_flat
  }.s2mPipe()

  io.skip_cmd << s2
}
