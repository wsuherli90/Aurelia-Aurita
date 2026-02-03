package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

class LaplacianOrder4Unit(config: AureliaConfig) extends Component {
  val io = new Bundle {
    val window_3d = in(Window5D()) 
    val g_inv     = in(Matrix3())
    
    val input_valid = in Bool()
    val output_ready = in Bool() 
    val output_laplacian = master Stream(Matrix3())
  }
  
  case class LapInput() extends Bundle {
     val w     = Window5D()
     val g_inv = Matrix3()
  }
  
  val input_stream = Stream(LapInput())
  input_stream.valid := io.input_valid
  input_stream.payload.w := io.window_3d
  input_stream.payload.g_inv := io.g_inv

  val INV_12H2 = RealConst(1.0 / (12.0 * config.gridDx * config.gridDx))
  val INV_4H2  = RealConst(1.0 / (4.0 * config.gridDx * config.gridDx))

  def calc_d2_pure_order4(p_plus2: Matrix3, p_plus1: Matrix3, p_center: Matrix3, p_minus1: Matrix3, p_minus2: Matrix3): Matrix3 = {
      val res = Matrix3()
      for(r <- 0 until 3; c <- 0 until 3) {
          val term = (-p_plus2(r)(c) + 
                      (RealConst(16.0) * p_plus1(r)(c)) - 
                      (RealConst(30.0) * p_center(r)(c)) + 
                      (RealConst(16.0) * p_minus1(r)(c)) - 
                      p_minus2(r)(c)) 
          res(r)(c) := term * INV_12H2
      }
      res
  }
  
  def calc_mixed_xy(w_plane: Vec[Vec[Matrix3]]): Matrix3 = {
     val res = Matrix3()
     val f_pp = w_plane(3)(3)
     val f_pm = w_plane(1)(3)
     val f_mp = w_plane(3)(1)
     val f_mm = w_plane(1)(1)
     
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
      
      p.dzz := calc_d2_pure_order4(w(4)(2)(2), w(3)(2)(2), w(2)(2)(2), w(1)(2)(2), w(0)(2)(2))
      
      p.dyy := calc_d2_pure_order4(w(2)(4)(2), w(2)(3)(2), w(2)(2)(2), w(2)(1)(2), w(2)(0)(2))
      
      p.dxx := calc_d2_pure_order4(w(2)(2)(4), w(2)(2)(3), w(2)(2)(2), w(2)(2)(1), w(2)(2)(0))
      
      p.dxy := calc_mixed_xy(w(2))
      
      val slice_yz = Vec(w.map(plane => plane.map(row => row(2)))) 
      
      p.dyz := calc_mixed_xy(slice_yz) 
      
      val slice_zx = Vec(w.map(plane => plane(2))) 
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
