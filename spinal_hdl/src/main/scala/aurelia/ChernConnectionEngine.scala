package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._

// Placeholder for NonLinear Connection
class NonLinearConnectionN extends Component {
  val io = new Bundle {
      // Connect as pure flow logic
      val gamma_osc = in(Tensor3())
      val y_vec     = in(Vector3())
      val N_matrix  = out(Matrix3())
  }

  // N^i_j = Sum_k (Gamma^i_jk * y^k)
  for(i <- 0 until 3) {
      for(j <- 0 until 3) {
          var sum = RealConst(0.0)
          for(k <- 0 until 3) {
              val term = io.gamma_osc(i)(j)(k) * io.y_vec(k)
              sum = sum + term
          }
          io.N_matrix(i)(j) := sum
      }
  }
}

class ChernConnectionEngine extends Component {
  val io = new Bundle {
    val win_x       = in(StencilWindow())
    val win_y       = in(StencilWindow())
    val win_z       = in(StencilWindow())
    val g_inv_in    = in(Matrix3())
    val y_vec       = in(Vector3())

    // Stream inputs wrapper
    val input_valid = in Bool()
    val output_ready = in Bool()
    
    val gamma_final = master Stream(Tensor3())
  }
  
  case class ChernInput() extends Bundle {
     val win_x    = StencilWindow()
     val win_y    = StencilWindow()
     val win_z    = StencilWindow()
     val g_inv    = Matrix3()
     val y_vec    = Vector3()
  }
  
  val input_stream = Stream(ChernInput())
  input_stream.valid := io.input_valid
  input_stream.payload.win_x := io.win_x
  input_stream.payload.win_y := io.win_y
  input_stream.payload.win_z := io.win_z
  input_stream.payload.g_inv := io.g_inv_in
  input_stream.payload.y_vec := io.y_vec

  val INV_12H = RealConst(1.0 / (12.0 * GRID_DX_C))

  // 4th Order 1st Derivative
  def calc_deriv_order4(w: StencilWindow): Matrix3 = {
     val res = Matrix3()
     for(r <- 0 until 3; c <- 0 until 3) {
         // (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / 12h
         val val_p2 = w(4)(r)(c)
         val val_p1 = w(3)(r)(c)
         val val_m1 = w(1)(r)(c)
         val val_m2 = w(0)(r)(c)
         
         val d_val = (-val_p2 + (RealConst(8.0)*val_p1) - (RealConst(8.0)*val_m1) + val_m2) * INV_12H
         res(r)(c) := d_val
     }
     res
  }

  // Stage 1: Calculate Derivatives
  case class DerivPayload() extends Bundle {
      // We need derivatives for all components.
      // D_k g_ij where k = x,y,z
      // We store them as a Vec of Matrices? Or 3 separate matrices for logic clarity?
      val d_x_g = Matrix3()
      val d_y_g = Matrix3()
      val d_z_g = Matrix3()
      val g_inv = Matrix3()
      val y_vec = Vector3()
  }
  
  val s1 = input_stream.translateWith {
      val i = input_stream.payload
      val p = DerivPayload()
      p.d_x_g := calc_deriv_order4(i.win_x)
      p.d_y_g := calc_deriv_order4(i.win_y)
      p.d_z_g := calc_deriv_order4(i.win_z)
      p.g_inv := i.g_inv
      p.y_vec := i.y_vec
      p
  }.s2mPipe()

  // Stage 2: Calculate Osculating Gamma
  case class GammaPayload() extends Bundle {
      val gamma_osc = Tensor3()
      val y_vec     = Vector3()
  }
  
  val s2 = s1.translateWith {
      val p = s1.payload
      val gamma_osc = Tensor3()
      
      val derivs = Vec(p.d_x_g, p.d_y_g, p.d_z_g)
      
       // Loop over k, i, j
       // Gamma^k_ij = 0.5 * g^km * (d_j g_im + d_i g_jm - d_m g_ij)
      for(k <- 0 until 3; i <- 0 until 3; j <- 0 until 3) {
          var sum = RealConst(0.0)
          
          for(m <- 0 until 3) {
              val d_j_gim = derivs(j)(i)(m)
              val d_i_gjm = derivs(i)(j)(m)
              val d_m_gij = derivs(m)(i)(j)
              
              val term = d_j_gim + d_i_gjm - d_m_gij
              sum = sum + (p.g_inv(k)(m) * term)
          }
          gamma_osc(k)(i)(j) := sum * RealConst(0.5)
      }
      
      val res = GammaPayload()
      res.gamma_osc := gamma_osc
      res.y_vec     := p.y_vec
      res
  }.s2mPipe() // This is compute intensive (Tensor contraction), s2mPipe is crucial

  // Stage 3: Nonlinear Connection Addition
  // We can instantiate the module inside the map? 
  // Ideally, Modules are instantiated outside, or we just write the logic.
  // We will write the logic inline or use the helper class combinatorially.
  
  val s3 = s2.translateWith {
      val p = s2.payload
      
      // Instantiate logic
      val uNonLinear = new NonLinearConnectionN()
      uNonLinear.io.gamma_osc := p.gamma_osc
      uNonLinear.io.y_vec     := p.y_vec
      
      val gamma_final = Tensor3()
      for(k <- 0 until 3; i <- 0 until 3; j <- 0 until 3) {
           gamma_final(k)(i)(j) := p.gamma_osc(k)(i)(j) + uNonLinear.io.N_matrix(i)(j)
      }
      gamma_final
  }.s2mPipe()
  
  io.gamma_final << s3
}
