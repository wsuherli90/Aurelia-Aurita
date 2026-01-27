package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

/**
 * Laplacian Unit (Order 4)
 * Input: 3 Stencils (X, Y, Z axis) + Inverse Metric
 * Output: Laplacian Matrix
 */
class LaplacianOrder4Unit extends Component {
  val io = new Bundle {
    val window_3d = in(Window3D()) 
    val g_inv     = in(Matrix3())
    
    val input_valid = in Bool()
    val output_ready = in Bool() 
    val output_laplacian = master Stream(Matrix3())
  }
  
  case class LapInput() extends Bundle {
     val w     = Window3D()
     val g_inv = Matrix3()
  }
  
  val input_stream = Stream(LapInput())
  input_stream.valid := io.input_valid
  input_stream.payload.w := io.window_3d
  input_stream.payload.g_inv := io.g_inv

  val INV_4H2 = RealConst(1.0 / (4.0 * GRID_DX_C * GRID_DX_C))
  val INV_H2 = RealConst(1.0 / (GRID_DX_C * GRID_DX_C))

  def calc_d2_pure(p_minus: Matrix3, p_center: Matrix3, p_plus: Matrix3): Matrix3 = {
      val res = Matrix3()
      for(r <- 0 until 3; c <- 0 until 3) {
          res(r)(c) := (p_plus(r)(c) - (RealConst(2.0)*p_center(r)(c)) + p_minus(r)(c)) * INV_H2
      }
      res
  }
  
  def calc_mixed_xy(w_plane: Vec[Vec[Matrix3]]): Matrix3 = {
     val res = Matrix3()
     val f_pp = w_plane(0)(0) 
     val f_pm = w_plane(0)(2) 
     val f_mp = w_plane(2)(0) 
     val f_mm = w_plane(2)(2) 
     
     for(r <- 0 until 3; c <- 0 until 3) {
         res(r)(c) := (f_pp(r)(c) - f_pm(r)(c) - f_mp(r)(c) + f_mm(r)(c)) * INV_4H2
     }
     res
  }
  
  case class DerivPayload() extends Bundle {
      val dxx = Matrix3(); val dyy = Matrix3(); val dzz = Matrix3()
      val dxy = Matrix3(); val dyz = Matrix3(); val dzx = Matrix3()
      
      val g_inv = Matrix3()
  }
  
  val s1 = input_stream.translateWith {
      val w = input_stream.payload.w
      val p = DerivPayload()
      p.g_inv := input_stream.payload.g_inv
      
      p.dxx := calc_d2_pure(w(1)(1)(2), w(1)(1)(1), w(1)(1)(0))
      p.dyy := calc_d2_pure(w(1)(2)(1), w(1)(1)(1), w(1)(0)(1))
      p.dzz := calc_d2_pure(w(2)(1)(1), w(1)(1)(1), w(0)(1)(1))
      
      p.dxy := calc_mixed_xy(w(1)) 
      
      val slice_yz = Vec(w(0).map(_(1)), w(1).map(_(1)), w(2).map(_(1))) 
      p.dyz := calc_mixed_xy(slice_yz)
      
      val slice_zx = Vec(w(0)(1), w(1)(1), w(2)(1))
      p.dzx := calc_mixed_xy(slice_zx)
      
      p
  }.s2mPipe() 
  
  val s2 = s1.translateWith {
      val p = s1.payload
      val lap = Matrix3()
      val g = p.g_inv
      
      for(a <- 0 until 3; b <- 0 until 3) {
          val term_ii = (g(0)(0) * p.dxx(a)(b)) + 
                        (g(1)(1) * p.dyy(a)(b)) + 
                        (g(2)(2) * p.dzz(a)(b))
                        
          val term_xy = g(0)(1) * p.dxy(a)(b)
          val term_yz = g(1)(2) * p.dyz(a)(b)
          val term_zx = g(2)(0) * p.dzx(a)(b)
          
          lap(a)(b) := term_ii + (RealConst(2.0) * (term_xy + term_yz + term_zx))
      }
      lap
  }.s2mPipe()
  
  io.output_laplacian << s2
}
