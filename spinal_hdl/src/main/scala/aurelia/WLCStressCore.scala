package aurelia

import spinal.core._
import spinal.lib._
import AureliaTypes._


// WLC Force-Extension: f ~ 1/4(1-x)^2 - 1/4 + x
// x = extension / L_contour. 
// We assume 'y_vec' length represents extension (or scaled).
// We calculate lambda = sqrt( g_ab y^a y^b ).
// x = lambda (assuming normalized units).
// Singularity at x=1. Max extension should be clamped.

class WLCStressCore extends Component {
  val io = new Bundle {
    val input_valid = in Bool()
    val output_ready = in Bool()
    
    val mu_in     = in(RealType())
    val g_metric  = in(Matrix3())
    val y_vec     = in(Vector3()) 
    
    val output = master Stream(Matrix3())
  }
  
  // Define input stream explicitly to manage backpressure
  case class WLCInput() extends Bundle {
     val mu = RealType()
     val g  = Matrix3()
     val y  = Vector3()
  }
  
  val inputStream = Stream(WLCInput())
  inputStream.valid := io.input_valid
  inputStream.payload.mu := io.mu_in
  inputStream.payload.g  := io.g_metric
  inputStream.payload.y  := io.y_vec
  
  // Pipeline Step 1: Compute Metric norm squared (lambda^2)
  // lambda^2 = g_ab * y^a * y^b
  
  case class S1Payload() extends Bundle {
     val lambda_sq = RealType()
     val mu        = RealType()
     val g         = Matrix3()
     val y         = Vector3()
     val y_outer   = Matrix3()
  }

  val s1 = inputStream.translateWith {
      val i = inputStream.payload
      val res = S1Payload()
      res.mu := i.mu
      res.g  := i.g
      res.y  := i.y
      
      val y_cov = Vector3()
      for(r <- 0 until 3) {
          y_cov(r) := (i.g(r)(0) * i.y(0)) + (i.g(r)(1) * i.y(1)) + (i.g(r)(2) * i.y(2))
      }
      
      // lambda^2 = y_a * y^a_cov
      res.lambda_sq := (y_cov(0) * i.y(0)) + (y_cov(1) * i.y(1)) + (y_cov(2) * i.y(2))
      
      // Precompute outer product y_a * y_b for later stress tensor
      for(r <- 0 until 3; c <- 0 until 3) {
          res.y_outer(r)(c) := i.y(r) * i.y(c)
      }
      res
  }.s2mPipe()
  
  // Pipeline Step 2: Sqrt(lambda^2) -> lambda
  // We use our new MathUtils.FixedSqrtPipelined
  
  val sqrtUnit = new MathUtils.FixedSqrtPipelined()
  
  // Fork logic: we need to send lambda_sq to Sqrt, but keep other data context bypassing.
  // Standard StreamJoin pattern.
  
  val (to_sqrt, context_bypass) = StreamFork2(s1)
  
  sqrtUnit.io.input << to_sqrt.translateWith(to_sqrt.payload.lambda_sq)
  
  val joined_lambda = StreamJoin.arg(sqrtUnit.io.output, context_bypass)
  
  // Pipeline Step 3: WLC Factor Calculation
  // f(x) = 0.25/(1-x)^2 - 0.25 + x
  // We need 1/(1-x).
  
  case class S3Payload() extends Bundle {
     val one_minus_x = RealType()
     val x           = RealType()
     val ctx         = S1Payload()
  }
  
  val s3 = joined_lambda.translateWith { 
      val lambda = joined_lambda.payload._1
      val ctx    = joined_lambda.payload._2
      
      val res = S3Payload()
      res.x   := lambda
      res.ctx := ctx
      
      // Clamp x to 0.95 to avoid singularity explosion
      val x_clamped = (lambda > RealConst(0.95)) ? RealConst(0.95) | lambda
      res.one_minus_x := RealConst(1.0) - x_clamped
      res
  }.s2mPipe()
  
  // Pipeline Step 4: Reciprocal
  
  val recipUnit = new MathUtils.ReciprocalUnit() 
  val (to_recip, ctx_s3) = StreamFork2(s3)
  
  recipUnit.io.input << to_recip.translateWith(to_recip.payload.one_minus_x)
  
  val recip_out = recipUnit.io.output

  // We need to synchronize the context to match ReciprocalUnit latency.
  // Since ReciprocalUnit is now dynamic (Stream), fixed delay is invalid.
  // We use a Queue (FIFO) to buffer the context until the reciprocal result is ready.
  // Depth 16 covers the ~6-10 cycle latency of NR iterations comfortably.
  
  val ctx_buffered = ctx_s3.queue(16)
  
  val join_wlc = StreamJoin.arg(recip_out, ctx_buffered)
  
  // Pipeline Step 5: Final Stress Tensor
  // T_ab = mu * (g_ab + wlc_factor * y_a * y_b)
  // wlc_factor = 0.25 * (rec^2) - 0.25 + x
  
  val s5 = join_wlc.translateWith {
      val rec = join_wlc.payload._1
      val p   = join_wlc.payload._2 // S3Payload
      
      val rec_sq = rec * rec
      val term1  = rec_sq * RealConst(0.25)
      val wlc_f  = term1 - RealConst(0.25) + p.x
      
      val res_stress = Matrix3()
      
      // Isotropic part + Active Anisotropic part
      // T = mu * (g + wlc * y \otimes y)
      // Usually active stress is added to stress tensor.
      // Here we output the "Active Stress Tensor".
      
      for(i <- 0 until 3; j <- 0 until 3) {
          val anisotropy = wlc_f * p.ctx.y_outer(i)(j)
          val total      = p.ctx.g(i)(j) + anisotropy
          res_stress(i)(j) := p.ctx.mu * total
      }
      res_stress
  }.s2mPipe()
  
  io.output << s5
  
  // Handshake management updates
  // ReciprocalUnit: we need to throttle input if output queue is full?
  // Since ReciprocalUnit has no 'ready', we must rely on the fact that if downstream is ready,
  // we are fine. If downstream is NOT ready, StreamJoin holds.
  // BUT: `ctx_s3.delay(6)` will backpressure!
  // If `ctx_delayed` is not consumed (join_wlc stalled), `ctx_s3` stalls.
  // If `ctx_s3` stalls, `s3` stalls.
  // If `s3` stalls, we stop feeding `recipUnit`.
  // Wait: `to_recip` is forked. If `ctx_s3` stalls, `to_recip` also stalls because Fork awaits both?
  // Yes, StreamFork awaits all branches to fire.
  // So if we don't fire `to_recip`, we don't feed `ReciprocalUnit`.
  // This assumes `to_recip.ready` is wired to logic that says "Can I feed recip logic?".
  // Since `ReciprocalUnit` is always consuming if we toggle valid, we can set `to_recip.ready := True`.
  // BUT valid loop: 
  // We should only fire `s3` if we can effectively push to the implicit pipeline. 
  // It works fine. `StreamFork` default handles this synchronization.
}
