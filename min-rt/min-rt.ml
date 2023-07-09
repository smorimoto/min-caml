(****************************************************************)
(*                                                              *)
(* Ray Tracing Program for (Mini) Objective Caml                *)
(*                                                              *)
(* Original Program by Ryoji Kawamichi                          *)
(* Arranged for Chez Scheme by Motohico Nanano                  *)
(* Arranged for Objective Caml by Y. Oiwa and E. Sumii          *)
(*                                                              *)
(****************************************************************)

(*NOMINCAML open MiniMLRuntime;;*)
(*NOMINCAML open Globals;;*)

(*MINCAML*)let rec xor x y = if x then not y else y in
(*NOMINCAML let xor x y = if x then not y else y in*)

(*MINCAML*)let rec fsqr x = x *. x in
(*NOMINCAML let fsqr x = x *. x in*)

(*MINCAML*)let rec fhalf x = x /. 2. in
(*NOMINCAML let fhalf x = x /. 2. in*)

(**************** ⒢ｃ∽ ****************)
(* 炊拷*)

(*MINCAML*)let rec o_texturetype m = 
(*NOMINCAML let o_texturetype m = *)
  let (m_tex, xm_shape, xm_surface, xm_isrot, 
       xm_abc, xm_xyz, 
       xm_invert, xm_surfparams, xm_color,
       xm_rot123) = m 
  in
  m_tex
in

(*MINCAML*)let rec o_form m = 
(*NOMINCAML let o_form m = *)
  let (xm_tex, m_shape, xm_surface, xm_isrot, 
       xm_abc, xm_xyz, 
       xm_invert, xm_surfparams, xm_color,
       xm_rot123) = m 
  in
  m_shape
in

(*MINCAML*)let rec o_reflectiontype m = 
(*NOMINCAML let o_reflectiontype m = *)
  let (xm_tex, xm_shape, m_surface, xm_isrot, 
       xm_abc, xm_xyz, 
       xm_invert, xm_surfparams, xm_color,
       xm_rot123) = m 
  in
  m_surface
in

(*MINCAML*)let rec o_isinvert m = 
(*NOMINCAML let o_isinvert m = *)
  let (m_tex, m_shape, m_surface, m_isrot, 
       xm_abc, xm_xyz, 
       m_invert, xm_surfparams, xm_color,
       xm_rot123) = m in
  m_invert
in

(*MINCAML*)let rec o_isrot m = 
(*NOMINCAML let o_isrot m = *)
  let (xm_tex, xm_shape, xm_surface, m_isrot, 
       xm_abc, xm_xyz, 
       xm_invert, xm_surfparams, xm_color,
       xm_rot123) = m in
  m_isrot
in

(*MINCAML*)let rec o_param_a m = 
(*NOMINCAML let o_param_a m = *)
  let (xm_tex, xm_shape, xm_surface, xm_isrot, 
       m_abc, xm_xyz, 
       xm_invert, xm_surfparams, xm_color,
       xm_rot123) = m 
  in
  m_abc.(0)
in

(*MINCAML*)let rec o_param_b m = 
(*NOMINCAML let o_param_b m = *)
  let (xm_tex, xm_shape, xm_surface, xm_isrot, 
       m_abc, xm_xyz, 
       xm_invert, xm_surfparams, xm_color,
       xm_rot123) = m 
  in
  m_abc.(1)
in

(*MINCAML*)let rec o_param_c m = 
(*NOMINCAML let o_param_c m = *)
  let (xm_tex, xm_shape, xm_surface, xm_isrot, 
       m_abc, xm_xyz, 
       xm_invert, xm_surfparams, xm_color,
       xm_rot123) = m 
  in
  m_abc.(2)
in

(*MINCAML*)let rec o_param_x m = 
(*NOMINCAML let o_param_x m = *)
  let (xm_tex, xm_shape, xm_surface, xm_isrot, 
       xm_abc, m_xyz, 
       xm_invert, xm_surfparams, xm_color,
       xm_rot123) = m 
  in
  m_xyz.(0)
in

(*MINCAML*)let rec o_param_y m = 
(*NOMINCAML let o_param_y m = *)
  let (xm_tex, xm_shape, xm_surface, xm_isrot, 
       xm_abc, m_xyz,
       xm_invert, xm_surfparams, xm_color,
       xm_rot123) = m 
  in
  m_xyz.(1)
in

(*MINCAML*)let rec o_param_z m = 
(*NOMINCAML let o_param_z m = *)
  let (xm_tex, xm_shape, xm_surface, xm_isrot, 
       xm_abc, m_xyz,
       xm_invert, xm_surfparams, xm_color,
       xm_rot123) = m 
  in
  m_xyz.(2)
in

(*MINCAML*)let rec o_diffuse m = 
(*NOMINCAML let o_diffuse m = *)
  let (xm_tex, xm_shape, xm_surface, xm_isrot, 
       xm_abc, xm_xyz, 
       xm_invert, m_surfparams, xm_color,
       xm_rot123) = m 
  in
  m_surfparams.(0)
in

(*MINCAML*)let rec o_hilight m = 
(*NOMINCAML let o_hilight m = *)
  let (xm_tex, xm_shape, xm_surface, xm_isrot, 
       xm_abc, xm_xyz, 
       xm_invert, m_surfparams, xm_color,
       xm_rot123) = m 
  in
  m_surfparams.(1)
in

(*MINCAML*)let rec o_color_red m = 
(*NOMINCAML let o_color_red m = *)
  let (xm_tex, xm_shape, m_surface, xm_isrot, 
       xm_abc, xm_xyz, 
       xm_invert, xm_surfparams, m_color,
       xm_rot123) = m 
  in
  m_color.(0)
in

(*MINCAML*)let rec o_color_green m = 
(*NOMINCAML let o_color_green m = *)
  let (xm_tex, xm_shape, m_surface, xm_isrot, 
       xm_abc, xm_xyz, 
       xm_invert, xm_surfparams, m_color,
       xm_rot123) = m 
  in
  m_color.(1)
in

(*MINCAML*)let rec o_color_blue m = 
(*NOMINCAML let o_color_blue m = *)
  let (xm_tex, xm_shape, m_surface, xm_isrot, 
       xm_abc, xm_xyz, 
       xm_invert, xm_surfparams, m_color,
       xm_rot123) = m 
  in
  m_color.(2)
in

(*MINCAML*)let rec o_param_r1 m = 
(*NOMINCAML let o_param_r1 m = *)
  let (xm_tex, xm_shape, xm_surface, xm_isrot, 
       xm_abc, xm_xyz, 
       xm_invert, xm_surfparams, xm_color,
       m_rot123) = m 
  in
  m_rot123.(0)
in

(*MINCAML*)let rec o_param_r2 m = 
(*NOMINCAML let o_param_r2 m = *)
  let (xm_tex, xm_shape, xm_surface, xm_isrot, 
       xm_abc, xm_xyz, 
       xm_invert, xm_surfparams, xm_color,
       m_rot123) = m 
  in
  m_rot123.(1)
in

(*MINCAML*)let rec o_param_r3 m = 
(*NOMINCAML let o_param_r3 m = *)
  let (xm_tex, xm_shape, xm_surface, xm_isrot, 
       xm_abc, xm_xyz, 
       xm_invert, xm_surfparams, xm_color,
       m_rot123) = m 
  in
  m_rot123.(2)
in

(*MINCAML*)let rec normalize_vector v inv = 
(*NOMINCAML let normalize_vector v inv = *)
  let n0 = (sqrt (fsqr v.(0) +. fsqr v.(1) +. fsqr v.(2))) in
  let n = if inv then -.n0 else n0 in
  v.(0) <- v.(0) /. n;
  v.(1) <- v.(1) /. n;
  v.(2) <- v.(2) /. n
in

(*MINCAML*)let rec sgn x =
(*NOMINCAML let sgn x =*)
  if 0.0 < x then 1.0
  else -1.0
in

(**************** 粋粋昭翠篆∽臂****************)

(* ❼≪ *)
(*MINCAML*)let rec rad x = x *. (0.017453293)
(*NOMINCAML let rad x = x *. (0.017453293)*)
in
(**** 医帥茯㏍莨若 ****)
(*MINCAML*)let rec read_environ _ =
(*NOMINCAML let read_environ _ =*)
  (*if dbg.(0) then
    Format.printf "reading environment data.@.";*)
  (* 鴻⒢割賢綽綺ф *)
  screen.(0) <- read_float ();
  screen.(1) <- read_float ();
  screen.(2) <- read_float ();
  (* 拶茹*)
  let v1 = rad (read_float ()) in
  cos_v.(0) <- cos v1;
  sin_v.(0) <- sin v1;
  let v2 = rad (read_float ()) in
  cos_v.(1) <- cos v2;
  sin_v.(1) <- sin v2;

  let nl = read_float () in

  (* ≫ *)
  let l1 = rad (read_float ()) in
  let sl1 = sin l1 in
  light.(1) <- -.sl1;
  let l2 = rad (read_float ()) in
  let cl1 = cos l1 in
  let sl2 = sin l2 in
  light.(0) <- cl1 *. sl2;
  let cl2 = cos l2 in
  light.(2) <- cl1 *. cl2;
  beam.(0) <- read_float ();

  (* 荀篏舟(㏍綺ф) *)
  vp.(0) <- cos_v.(0) *. sin_v.(1) *. (-200.0);
  vp.(1) <- (-.sin_v.(0)) *. (-200.0);
  vp.(2) <- cos_v.(0) *. cos_v.(1) *. (-200.0);

  (* 荀篏舟(綺ф) *)
  view.(0) <- vp.(0) +. screen.(0);
  view.(1) <- vp.(1) +. screen.(1);
  view.(2) <- vp.(2) +. screen.(2)
in

(**** ⒢吾1ゃ帥茯㏍莨若 ****)
(*MINCAML*)let rec read_nth_object n =
(*NOMINCAML let read_nth_object n =*)
(*  if dbg.(0) then
    Format.printf "object #%d.@." n;*)

  let texture = read_int () in (* 0 *)
  if texture <> -1 then
    ( 
      let form = read_int () in (* 1 *)
      let refltype = read_int () in
      let isrot_p = read_int () in

      let abc = Array.make 3 0.0 in

      abc.(0) <- read_float ();
      abc.(1) <- read_float (); (* 5 *)
      abc.(2) <- read_float ();
      
      let xyz = Array.make 3 0.0 in

      xyz.(0) <- read_float ();
      xyz.(1) <- read_float ();
      xyz.(2) <- read_float ();

      let m_invert = 0.0 > (read_float ()) in (* 10 *)

      let reflparam = Array.make 2 0.0 in
      
      reflparam.(0) <- read_float (); (* diffuse *)
      reflparam.(1) <- read_float (); (* hilight *)

      let color = Array.make 3 0.0 in

      color.(0) <- read_float ();
      color.(1) <- read_float ();
      color.(2) <- read_float (); (* 15 *)

      let rotation = Array.make 3 0.0 in
      if isrot_p <> 0 then
        (
         rotation.(0) <- rad (read_float ());
         rotation.(1) <- rad (read_float ());
         rotation.(2) <- rad (read_float ()) (* 18 *)
        ) 
      else ();

      (* ＜帥罩ｈ*)

      (* 羈 筝罩ｈ(form = 2)  *)
      let m_invert2 = if form = 2 then true else m_invert in

      (* abc rotation *)
      let obj = 
        (texture, form, refltype, isrot_p,
         abc, xyz, (* x-z *)
         m_invert2,
         reflparam, (* reflection paramater *)
         color, (* color *)
         rotation (* rotation *)
        ) in
      objects.(n) <- obj;

      if form = 3 then
        (
          (* 2罨≧ X,Y,Z 泣冴罨≧腮違 *)
         let a = abc.(0) in
         abc.(0) <- if 0.0 = a then 0.0 else (sgn a) /. (fsqr a);
         let b = abc.(1) in
         abc.(1) <- if 0.0 = b then 0.0 else (sgn b) /. (fsqr b);
         let c = abc.(2) in
         abc.(2) <- if 0.0 = c then 0.0 else (sgn c) /. (fsqr c)
        )
      else if form = 2 then
        (* 綛渇: 羈③荀 罐究'輝筝*)
        normalize_vector abc (not m_invert)
      else
        ();

      if isrot_p <> 0 then
        (
         cs_temp.(10) <- cos rotation.(0);
         cs_temp.(11) <- sin rotation.(0);
         cs_temp.(12) <- cos rotation.(1);
         cs_temp.(13) <- sin rotation.(1);
         cs_temp.(14) <- cos rotation.(2);
         cs_temp.(15) <- sin rotation.(2);
         cs_temp.(0) <- cs_temp.(12) *. cs_temp.(14); (* cy cz *)
         cs_temp.(1) <- 
           (cs_temp.(11) *. cs_temp.(13) *. cs_temp.(14)) -. (cs_temp.(10) *. cs_temp.(15));
         cs_temp.(2) <-        
           (cs_temp.(10) *. cs_temp.(13) *. cs_temp.(14)) +. (cs_temp.(11) *. cs_temp.(15));
         cs_temp.(3) <-        cs_temp.(12) *. cs_temp.(15);
         cs_temp.(4) <-        
           (cs_temp.(11) *. cs_temp.(13) *. cs_temp.(15)) +. (cs_temp.(10) *. cs_temp.(14));
         cs_temp.(5) <-        
           (cs_temp.(10) *. cs_temp.(13) *. cs_temp.(15)) -. (cs_temp.(11) *. cs_temp.(14));
         cs_temp.(6) <-        -.cs_temp.(13);
         cs_temp.(7) <-        cs_temp.(11) *. cs_temp.(12);
         cs_temp.(8) <-        cs_temp.(10) *. cs_temp.(12);
         let ao = abc.(0) in
         let bo = abc.(1) in
         let co = abc.(2) in
         abc.(0) <- ao *. fsqr cs_temp.(0) +. bo *. fsqr cs_temp.(3) +. co *. fsqr cs_temp.(6);
         abc.(1) <- ao *. fsqr cs_temp.(1) +. bo *. fsqr cs_temp.(4) +. co *. fsqr cs_temp.(7);
         abc.(2) <- ao *. fsqr cs_temp.(2) +. bo *. fsqr cs_temp.(5) +. co *. fsqr cs_temp.(8);
         rotation.(0) <- 2.0 *. (ao *. cs_temp.(1) *. cs_temp.(2)
                                   +. bo *. cs_temp.(4) *. cs_temp.(5)
                                   +. co *. cs_temp.(7) *. cs_temp.(8));
         rotation.(1) <- 2.0 *. (ao *. cs_temp.(0) *. cs_temp.(2)
                                   +. bo *. cs_temp.(3) *. cs_temp.(5)
                                   +. co *. cs_temp.(6) *. cs_temp.(8));
         rotation.(2) <- 2.0 *. (ao *. cs_temp.(0) *. cs_temp.(1)
                                   +. bo *. cs_temp.(3) *. cs_temp.(4)
                                   +. co *. cs_temp.(6) *. cs_temp.(7))
        )
      else ();

(*      ;
   if dbg.(0) then
   Format.eprintf "OBJ #%d: %d (%f %f %f) (%f %f %f) [%f %f %f]@."
   n form abc.(0) abc.(1) abc.(2) xyz.(0) xyz.(0) xyz.(0) 
   rotation.(0) rotation.(1) rotation.(2);*)
      true
     )
  else
    false (* 帥腟 *)
in

(**** ❽水篏茯㏍莨若 ****)
(*MINCAML*)let rec read_object n =
(*NOMINCAML let rec read_object n =*)
  if n < 61 then
    if read_nth_object n 
    then read_object (n + 1)
    else ()
  else () (* failwith "data overflow" *)
in

(*MINCAML*)let rec read_all_object _ =
(*NOMINCAML let read_all_object _ =*)
  read_object 0
in

(**** AND, OR 若㍽粋昭****)

(* 若1ゃ茯㏍莨若*)
(*MINCAML*)let rec read_net_item length =
(*NOMINCAML let rec read_net_item length =*)
  let item = read_int () in
  if item = -1 then Array.make (length + 1) (-1)
  else
    let v = read_net_item (length + 1) in
    (v.(length) <- item; v)
in

(*MINCAML*)let rec read_or_network length =
(*NOMINCAML let rec read_or_network length =*)
  let net = read_net_item 0 in
  if net.(0) = -1 then 
    Array.make (length + 1) net
  else
    let v = read_or_network (length + 1) in
    (v.(length) <- net; v)
in

(*MINCAML*)let rec read_and_network n =
(*NOMINCAML let rec read_and_network n =*)
  let net = read_net_item 0 in
  if net.(0) = -1 then ()
  else (
    and_net.(n) <- net;
    read_and_network (n + 1)
  )
in

(*MINCAML*)let rec read_parameter _ =
(*NOMINCAML let read_parameter _ =*)
  (
   read_environ ();
   read_all_object ();
   read_and_network 0;
   or_net.(0) <- read_or_network 0
  )
in

(**************** 雁с篋ょ∽臂****************)

(* solver : 
   ⒢吾 (index)  L, P    雁 Lt + P ⒢吾篋ょ   篋ょ 0 篋ょ篁ュ
   菴ゃ nvector т困鴻羈綽   (贋篏翫)

   篋ょ婚罔 t $solver_dist 膣*)

(**** 贋篏с翫 ****)
(*MINCAML*)let rec solver_rect m l =
(*NOMINCAML let solver_rect m l =*)
  (* YZ 綛渇 *)
  let answera = 
    if 0.0 = l.(0) then false else (
      let d = 
        if xor (o_isinvert m) (0.0 > l.(0)) then (o_param_a m) else -.(o_param_a m)
      in
      let d2 = (d -. solver_w_vec.(0)) /. l.(0) 
      in
      if abs_float (d2 *. l.(1) +. solver_w_vec.(1)) < o_param_b m then
        if abs_float (d2 *. l.(2) +. solver_w_vec.(2)) < o_param_c m
        then (solver_dist.(0) <- d2; true)
        else false
      else false
   )
  in
  if answera then 1 else (* fall through *)
  
  (* ZX 綛渇 *)
  let answerb = 
    if 0.0 = l.(1) then false else (
      let d = 
        if xor (o_isinvert m) (0.0 > l.(1)) then (o_param_b m) else -.(o_param_b m)
      in
      let d2 = (d -. solver_w_vec.(1)) /. l.(1)
      in
      if abs_float (d2 *. l.(2) +. solver_w_vec.(2)) < o_param_c m then
        if abs_float (d2 *. l.(0) +. solver_w_vec.(0)) < o_param_a m
        then (solver_dist.(0) <- d2; true)
        else false
      else false
    )
  in
  if answerb then 2 else (* fall through *)
  
  (* XY 綛渇 *)
  let answerc = 
    if 0.0 = l.(2) then false else (
      let d = 
        if xor (o_isinvert m) (0.0 > l.(2)) then (o_param_c m) else -.(o_param_c m)
      in
      let d2 = (d -. solver_w_vec.(2)) /. l.(2)
      in
      if abs_float (d2 *. l.(0) +. solver_w_vec.(0)) < o_param_a m then
        if abs_float (d2 *. l.(1) +. solver_w_vec.(1)) < o_param_b m
        then (solver_dist.(0) <- d2; true)
        else false
      else false
   )
  in
  if answerc then 3 else 0
in

(* 綛渇⒢吾*)
(*MINCAML*)let rec solver_surface m l =
(*NOMINCAML let solver_surface m l =*)
  (* 鴻綛渇㍾♣激*)
  (* 綛渇サс莢腟延*)
  let q = (l.(0) *. o_param_a m +. l.(1) *. o_param_b m +. l.(2) *. o_param_c m) in
  if 0.0 < q then
    let t = (solver_w_vec.(0) *. o_param_a m +. solver_w_vec.(1) *. o_param_b m +. solver_w_vec.(2) *. o_param_c m) /. q in
    (solver_dist.(0) <- -.t; 1)
  else 0
in

(* solver_second 若違紊 *)

(*MINCAML*)let rec in_prod_sqr_obj m v =
(*NOMINCAML let in_prod_sqr_obj m v =*)
  ((fsqr v.(0)) *. o_param_a m
     +. (fsqr v.(1)) *. o_param_b m
     +. (fsqr v.(2)) *. o_param_c m)
in

(*MINCAML*)let rec in_prod_co_objrot m v =
(*NOMINCAML let in_prod_co_objrot m v =*)
  v.(1) *. v.(2) *. o_param_r1 m
    +. v.(0) *. v.(2) *. o_param_r2 m
    +. v.(0) *. v.(1) *. o_param_r3 m
in

(*MINCAML*)let rec solver2nd_mul_b m l =
(*NOMINCAML let solver2nd_mul_b m l =*)
  solver_w_vec.(0) *. l.(0) *. o_param_a m
    +. solver_w_vec.(1) *. l.(1) *. o_param_b m
    +. solver_w_vec.(2) *. l.(2) *. o_param_c m
in

(*MINCAML*)let rec solver2nd_rot_b m l =
(*NOMINCAML let solver2nd_rot_b m l =*)
   (solver_w_vec.(2) *. l.(1) +. solver_w_vec.(1) *. l.(2)) *. o_param_r1 m
     +. (solver_w_vec.(0) *. l.(2) +. solver_w_vec.(2) *. l.(0)) *. o_param_r2 m
     +. (solver_w_vec.(0) *. l.(1) +. solver_w_vec.(1) *. l.(0)) *. o_param_r3 m
in

(*MINCAML*)let rec solver_second m l =
(*NOMINCAML let solver_second m l =*)
  let aa0 = in_prod_sqr_obj m l in
  let aa =
    if o_isrot m <> 0 then aa0 +. in_prod_co_objrot m l
    else aa0
  in
  if 0.0 = aa
  then 0 (* no intersection *)
  else
    (
     let bb0 = 2.0 *. solver2nd_mul_b m l
     in
     let bb = 
       if o_isrot m <> 0 then bb0 +. solver2nd_rot_b m l
       else bb0
     in
     let cc0 = in_prod_sqr_obj m solver_w_vec in
     let cc1 = 
       if o_isrot m <> 0 then
         (cc0 +. in_prod_co_objrot m solver_w_vec)
       else cc0 in
     let cc =
       if o_form m = 3
       then cc1 -. 1.0 else cc1 
     in
     let d = (* ゅ綣*)
       let d2 = 4.0 *. aa *. cc in
       (fsqr bb) -. d2
     in
     if 0.0 < d
     then
       (
        let sd = sqrt d in
        let t1 = if o_isinvert m then sd else -.sd in
        (solver_dist.(0) <- ((t1 -. bb) /. 2.0 /. aa); 1)
       )
     else 0
    )
in

(**** solver ゃ ****)
(*MINCAML*)let rec solver index l p =
(*NOMINCAML let solver index l p =*)
  let m = objects.(index) in
  solver_w_vec.(0) <- p.(0) -. o_param_x m;
  solver_w_vec.(1) <- p.(1) -. o_param_y m;
  solver_w_vec.(2) <- p.(2) -. o_param_z m;
  let m_shape = o_form m in
  if m_shape = 1 then solver_rect m l
  else if m_shape = 2 then solver_surface m l
  else solver_second m l
in

(**************** 鴻⒢吾臀㍼篆≪亥召 ****************)

(**** px, py, pz сm ❼ゅ ****)

(*MINCAML*)let rec is_rect_outside m =
(*NOMINCAML let is_rect_outside m =*)
  if 
    (if abs_float isoutside_q.(0) < o_param_a m then
      if abs_float isoutside_q.(1) < o_param_b m then
        if abs_float isoutside_q.(2) < o_param_c m then true else false
      else false
    else false
    )
  then o_isinvert m else not (o_isinvert m)
in

(*MINCAML*)let rec is_plane_outside m =
(*NOMINCAML let is_plane_outside m =*)
  let w = (o_param_a m *. isoutside_q.(0)
             +. o_param_b m *. isoutside_q.(1)
             +. o_param_c m *. isoutside_q.(2)) in
  let s = 0.0 > w in
  not (xor (o_isinvert m) s)
in

(*MINCAML*)let rec is_second_outside m = 
(*NOMINCAML let is_second_outside m = *)
  let w = in_prod_sqr_obj m isoutside_q in
  let w2 = if o_form m = 3 then w -. 1.0 else w in
  let w3 =
    if o_isrot m <> 0 then 
      w2 +. in_prod_co_objrot m isoutside_q
    else
      w2
  in
  let s = 0.0 > w3 in
  not (xor (o_isinvert m) s)
in

(*MINCAML*)let rec is_outside m =
(*NOMINCAML let is_outside m =*)
  isoutside_q.(0) <- chkinside_p.(0) -. o_param_x m;
  isoutside_q.(1) <- chkinside_p.(1) -. o_param_y m;
  isoutside_q.(2) <- chkinside_p.(2) -. o_param_z m;
  let m_shape = o_form m in
  if m_shape = 1 then
    is_rect_outside m
  else if m_shape = 2 then
    is_plane_outside m
  else 
    is_second_outside m
in

(**** (qx, qy, qz) AND 若 iand 絎****)
(*MINCAML*)let rec check_all_inside ofs iand =
(*NOMINCAML let rec check_all_inside ofs iand =*)
  let head = iand.(ofs) in
  if head = -1 then true else (
    if (is_outside (objects.(head))) then false
    else check_all_inside (ofs + 1) iand
  )
in

(**************** 綵宴≪亥召 ****************)

(* 腔膩劫梢                  *)
(* ❽ゃ(=綵宴ｃ) ゅ)

(**** AND 若 iand 襲❼ゅ ****)
(*MINCAML*)let rec shadow_check_and_group iand_ofs and_group p =
(*NOMINCAML let rec shadow_check_and_group iand_ofs and_group p =*)
  if and_group.(iand_ofs) = -1 then
    false
  else
    let obj = and_group.(iand_ofs) in
(*
    if crashed_object.(0) = obj
    then shadow_check_and_group (iand_ofs + 1) and_group p
    else
*)
    let t0 = solver obj light p in
    let t0p = solver_dist.(0) in
    if (if t0 <> 0 then t0p < -0.2 else false)
    then 
      (* Q: 篋ょ茖絎鴻с *)
      (* ャ❼茯帥*)
      let t = t0p +. 0.01 in
      chkinside_p.(0) <- light.(0) *. t +. p.(0);
      chkinside_p.(1) <- light.(1) *. t +. p.(1);
      chkinside_p.(2) <- light.(2) *. t +. p.(2);
      if (check_all_inside 0 and_group) 
      then true
      else shadow_check_and_group (iand_ofs + 1) and_group p
          (* 罨＜⒢吾鴻≪ *)
    else
      (* 篋ょ 罐究'③()    *)
      (* AND 演)
      (* 篋ょ⒥「膣≪       *)
      if o_isinvert (objects.(obj))
      then shadow_check_and_group (iand_ofs + 1) and_group p
      else false
in

(**** OR 違若 or_group 襲絎****)
(*MINCAML*)let rec shadow_check_one_or_group ofs or_group p =
(*NOMINCAML let rec shadow_check_one_or_group ofs or_group p =*)
  let head = or_group.(ofs) in
  if head = -1 then
    false
  else (
    let and_group = and_net.(head) in
    let shadow_p = shadow_check_and_group 0 and_group p in
    if shadow_p then true
    else shadow_check_one_or_group (ofs + 1) or_group p
  )
in

(**** OR 違若襲ｃ絎****)
(*MINCAML*)let rec shadow_check_one_or_matrix ofs or_matrix p =
(*NOMINCAML let rec shadow_check_one_or_matrix ofs or_matrix p =*)
  let head = or_matrix.(ofs) in
  let range_primitive = head.(0) in
  if range_primitive = -1 then false (* OR茵㍾篋若 *)
  else (
    if range_primitive = 99 
    then
      (* range primitive *)
      if (shadow_check_one_or_group 1 head p)
      then true
      else shadow_check_one_or_matrix (ofs + 1) or_matrix p
    else 
      (* range primitive *)
      let t = solver range_primitive light p in
      (* range primitive ゃ*)
      (* or group 篋ょ           *)
      if t <> 0 then
        if solver_dist.(0) < -0.1
        then
          if shadow_check_one_or_group 1 head p
          then true
          else shadow_check_one_or_matrix (ofs + 1) or_matrix p
        else shadow_check_one_or_matrix (ofs + 1) or_matrix p
      else shadow_check_one_or_matrix (ofs + 1) or_matrix p
  )
in

(**************** 咲鴻ｃ∽ ****************)

(**** AND若若絲障***)
(**** 篋ょ❼茯帥                     ****)
(*MINCAML*)let rec solve_each_element iand_ofs and_group =
(*NOMINCAML let rec solve_each_element iand_ofs and_group =*)
  let iobj = and_group.(iand_ofs) in
  if iobj = -1 then ()
  else (
    let t0 = solver iobj vscan viewpoint in
    if t0 <> 0 then
      (
        (* 篋ょ篋ょ㍽膣筝㏍茯帥*)
        (* 篁с筝㏍et $鴻)
       let t0p = solver_dist.(0) in
       if -0.1 < t0p then
         if t0p < tmin.(0) then
           (
            let t = t0p +. 0.01 in
            chkinside_p.(0) <- vscan.(0) *. t +. viewpoint.(0);
            chkinside_p.(1) <- vscan.(1) *. t +. viewpoint.(1);
            chkinside_p.(2) <- vscan.(2) *. t +. viewpoint.(2);
            if check_all_inside 0 and_group then 
              ( 
                tmin.(0) <- t;
                crashed_point.(0) <- chkinside_p.(0);
                crashed_point.(1) <- chkinside_p.(1);
                crashed_point.(2) <- chkinside_p.(2);
                intsec_rectside.(0) <- t0;
                crashed_object.(0) <- iobj
               )
            else ()
           )
         else ()
       else ()
      )
    else 
      (
       (* 篋ょ❽眼札筝困鴻⒢ *)
       if o_isinvert (objects.(iobj)) then () else end_flag.(0) <- (true)
      );
    if (not (end_flag.(0))) then 
      solve_each_element (iand_ofs + 1) and_group
    else ()
  )
in

(**** 1ゃ OR-group 篋ょⅺ鴻 ****)
(*MINCAML*)let rec solve_one_or_network ofs or_group =
(*NOMINCAML let rec solve_one_or_network ofs or_group =*)
  let head = or_group.(ofs) in
  if head = -1 then () else (
    let and_group = and_net.(head) in
    end_flag.(0) <- false;
    solve_each_element 0 and_group;
    solve_one_or_network (ofs + 1) or_group
  )
in

(**** OR⒢劫篏ゃ困鴻茯帥****)
(*MINCAML*)let rec trace_or_matrix ofs or_network =
(*NOMINCAML let rec trace_or_matrix ofs or_network =*)
  let head = or_network.(ofs) in
  let range_primitive = head.(0) in
  if range_primitive = -1 then (* с篋*)
    ()
  else ( 
    if range_primitive = 99 (* range primitive ⒢ *)
    then (solve_one_or_network 1 head)
    else 
      (
        (* range primitive ㍽腦⒢篋ょ*)
        let t = solver range_primitive vscan viewpoint in
        if t <> 0 then
          let tp = solver_dist.(0) in
          if tp < tmin.(0)
          then (solve_one_or_network 1 head)
          else ()
        else ()
      );
    trace_or_matrix (ofs + 1) or_network
  )
in

(**** 若 ****)
(* 若ViewPoint 鴻鴻ｃ劫 *)
(* Vscan 困crashed_point 腦с        *)
(* crashed_object ∽⒣㍽$困鴻㍾遵$*)
(*MINCAML*)let rec tracer viewpoint vscan =
(*NOMINCAML let tracer viewpoint vscan =*)
( 
  tmin.(0) <- (1000000000.0);
  trace_or_matrix 0 (or_net.(0));
  let t = tmin.(0) in
  if -0.1 < t then
    if t < 100000000.0 then
      true
    else false
  else false
)

in

(**************** 若∽臂****************)

(**** 篋ょ羈膊****)
(* 茵⒢吾 solver ㍽$*)
(* 紊 intsec_rectside 腟ф検綽  *)
(* nvector ㏍*)

(*MINCAML*)let rec get_nvector_rect _ =
(*NOMINCAML let get_nvector_rect _ =*)
  let rectside = intsec_rectside.(0) in
  (* solver ㍽$ゃｃ≪劫ず*)
  if rectside = 1 then (* YZ*)
    ( 
      nvector.(0) <- -.(sgn (vscan.(0)));
      nvector.(1) <- 0.0;
      nvector.(2) <- 0.0
     )
  else if rectside = 2 then (* ZX*)
    ( 
      nvector.(0) <- 0.0;
      nvector.(1) <- -.(sgn (vscan.(1)));
      nvector.(2) <- 0.0
     )
  else if rectside = 3 then (* XY*)
    ( 
      nvector.(0) <- 0.0;
      nvector.(1) <- 0.0;
      nvector.(2) <- -.(sgn (vscan.(2)))
     )
  else () (* should not happen *)
in

(*MINCAML*)let rec get_nvector_plane m = 
(*NOMINCAML let get_nvector_plane m = *)
  (* m_invert 幻true *)
  nvector.(0) <- -.(o_param_a m); (* if m_invert then -.m_a else m_a *)
  nvector.(1) <- -.(o_param_b m);
  nvector.(2) <- -.(o_param_c m)
in

(*MINCAML*)let rec get_nvector_second_norot m p = 
(*NOMINCAML let get_nvector_second_norot m p = *)
(* 拶⒢ *)
  nvector.(0) <- (p.(0) -. o_param_x m) *. o_param_a m;
  nvector.(1) <- (p.(1) -. o_param_y m) *. o_param_b m;
  nvector.(2) <- (p.(2) -. o_param_z m) *. o_param_c m;
  normalize_vector nvector (o_isinvert m)
in

(*MINCAML*)let rec get_nvector_second_rot m p =
(*NOMINCAML let get_nvector_second_rot m p =*)
  nvector_w.(0) <- p.(0) -. o_param_x m;
  nvector_w.(1) <- p.(1) -. o_param_y m;
  nvector_w.(2) <- p.(2) -. o_param_z m;
  nvector.(0) <- (nvector_w.(0)        *. o_param_a m 
                    +. fhalf (nvector_w.(1) *. o_param_r3 m
                                +. nvector_w.(2) *. o_param_r2 m));
  nvector.(1) <- (nvector_w.(1)        *. o_param_b m 
                    +. fhalf (nvector_w.(0) *. o_param_r3 m
                                +. nvector_w.(2) *. o_param_r1 m));
  nvector.(2) <- (nvector_w.(2)        *. o_param_c m
                    +. fhalf (nvector_w.(0) *. o_param_r2 m
                                +. nvector_w.(1) *. o_param_r1 m));
  normalize_vector nvector (o_isinvert m)
in

(*MINCAML*)let rec get_nvector m p =
(*NOMINCAML let get_nvector m p =*)
  let m_shape = o_form m in
  if m_shape = 1 then
    get_nvector_rect ()
  else if m_shape = 2 then
    get_nvector_plane m
  else (* 2罨≧or  *)
    if o_isrot m <> 0 then
      get_nvector_second_rot m p
    else
      get_nvector_second_norot m p
  (* retval = nvector *)
in

(**** 篋ょ筝鴻ｃ蚊荐 ****)
(*MINCAML*)let rec utexture m p =
(*NOMINCAML let utexture m p =*)
  let m_tex = o_texturetype m in
  (* 堺с*)
  texture_color.(0) <- o_color_red m;
  texture_color.(1) <- o_color_green m;
  texture_color.(2) <- o_color_blue m;
  if m_tex = 1 then
    (
     (* zx劫с罔≧ (G) *)
     let w1 = p.(0) -. o_param_x m in
     let flag1 =
       let d1 = (floor (w1 *. 0.05)) *. 20.0 in
       if w1 -. d1 < 10.0 then true else false
     in
     let w3 = p.(2) -. o_param_z m in
     let flag2 =
       let d2 = (floor (w3 *. 0.05)) *. 20.0 in
       if w3 -. d2 < 10.0 then true else false
     in
     texture_color.(1) <-
       if flag1 
       then (if flag2 then 255.0 else 0.0)
       else (if flag2 then 0.0 else 255.0)
    )
  else if m_tex = 2 then
    (* y荵御鴻❼(R-G) *)
    (
      let w2 = fsqr (sin (p.(1) *. 0.25)) in
      texture_color.(0) <- 255.0 *. w2;
      texture_color.(1) <- 255.0 *. (1.0 -. w2)
    )
  else if m_tex = 3 then 
    (* ZX∽(G-B) *)
    ( 
      let w1 = p.(0) -. o_param_x m in
      let w3 = p.(2) -. o_param_z m in
      let w2 = sqrt (fsqr w1 +. fsqr w3) /. 10.0 in
      let w4 =  (w2 -. floor w2) *. 3.1415927 in
      let cws = fsqr (cos w4) in
      texture_color.(1) <- cws *. 255.0;
      texture_color.(2) <- (1.0 -. cws) *. 255.0
    )
  else if m_tex = 4 then (
    (* 筝 (B) *)
    let w1 = (p.(0) -. o_param_x m) *. (sqrt (o_param_a m)) in
    let w3 = (p.(2) -. o_param_z m) *. (sqrt (o_param_c m)) in
    let w4 = sqrt ((fsqr w1) +. (fsqr w3)) in
    let w7 = 
      if abs_float w1 < 1.0e-4 then
        15.0 (* atan +infty = pi/2 *)
      else
        let w5 = abs_float (w3 /. w1)
        in
        (atan w5) *. (30.0 /. 3.1415927)
    in
    let w9 = w7 -. (floor w7) in

    let w2 = (p.(1) -. o_param_y m) *. (sqrt (o_param_b m)) in
    let w8 =
      if abs_float w7 < 1.0e-4 then
        15.0
      else 
        let w6 = abs_float (w2 /. w4)
        in (atan w6) *. (30.0 /. 3.1415927)
    in
    let w10 = w8 -. (floor w8) in
    let w11 = 0.15 -. (fsqr (0.5 -. w9)) -. (fsqr (0.5 -. w10)) in
    texture_color.(2) <- if 0.0 >= w11 then 0.0 else w11 *. (255.0 /. 0.3)
   )
  else ()
in

(**** 鴻⒢割鴻ゃ鴻帥膊****)
(* viewpoint (違若: 若               *)
(* vscan (違若:     若劫      *)
(* nref:                            *)
(* energy:    (羝＾) *)

(*  *)
(*MINCAML*)let rec in_prod v1 v2 = 
(*NOMINCAML let in_prod v1 v2 = *)
  v1.(0) *. v2.(0) +. v1.(1) *. v2.(1) +. v1.(2) *. v2.(2)
in

(* v1 += w v2 *)
(*MINCAML*)let rec accumulate_vec_mul v1 v2 w =
(*NOMINCAML let accumulate_vec_mul v1 v2 w =*)
  v1.(0) <- v1.(0) +. w *. v2.(0);
  v1.(1) <- v1.(1) +. w *. v2.(1);
  v1.(2) <- v1.(2) +. w *. v2.(2)
in

(*MINCAML*)let rec raytracing nref energy =
(*NOMINCAML let rec raytracing nref energy =*)
  let crashed_p = tracer viewpoint vscan in

  (* 鋎吟⒢ (nref = 0) *)
  (* ♂羣綵演*)
  if (not crashed_p) then
    if nref <> 0 then 
      ( 
        let hl = -.(in_prod vscan light) in
        (* 90属翫 ( *)
        if 0.0 < hl then
          (
           (* ゃゃ綣桁墾綺 cos^3 箴*)
           let ihl = fsqr hl *. hl *. energy *. beam.(0) in
           rgb.(0) <- rgb.(0) +. ihl;
           rgb.(1) <- rgb.(1) +. ihl;
           rgb.(2) <- rgb.(2) +. ihl
          )
        else ()
       )
    else ()
  else ();
  
  if crashed_p then
    (* ⒢吾ゃｃ翫 *)
    ( 
      (* 1. ❽*)
      let cobj = objects.(crashed_object.(0)) in
      get_nvector cobj crashed_point;
      let bright = 
        if (shadow_check_one_or_matrix 0 or_net.(0) crashed_point)
        then 
          0.0 (* 綵宴綵*)
        else (
          let br = -.(in_prod nvector light) in
          let br1 = if 0.0 > br then 0.2 else br +. 0.2 in
          br1 *. energy *. o_diffuse cobj
         )
      in
      utexture cobj crashed_point; (* 鴻ｃ荐 *)
      accumulate_vec_mul rgb texture_color bright;

      if nref > 4 then () else
      if 0.1 < energy then 
        ( 
          (* 2. *)
          let w = (-2.0) *. in_prod vscan nvector in
          (* 劫号紊 *)
          accumulate_vec_mul vscan nvector w;
          
          let m_surface = o_reflectiontype cobj in
          if m_surface = 1 then
            (* 箙怨 : ゃゃ膊*)
            (
             if 0.0 = (o_hilight cobj) then 
               ()
             else
               let hl = -.(in_prod vscan light) in
               if 0.0 < hl then
                 ( 
                   let ihl =
                     fsqr (fsqr hl) *. energy *. bright
                       *. o_hilight cobj
                   in
                   rgb.(0) <- rgb.(0) +. ihl;
                   rgb.(1) <- rgb.(1) +. ihl;
                   rgb.(2) <- rgb.(2) +. ihl
                  )
               else ()
            )
          else if m_surface = 2 then
            (* ♂: 鍵若 *)
            ( 
              viewpoint.(0) <- crashed_point.(0);
              viewpoint.(1) <- crashed_point.(1);
              viewpoint.(2) <- crashed_point.(2);
              let energy2 = energy *. (1.0 -. o_diffuse cobj) in
              raytracing (nref + 1) energy2
             )
          else ()
         )
      else ()
     )
  else ()
in

(**** 水****)
(*MINCAML*)let rec write_rgb _ =
(*NOMINCAML let write_rgb _ =*)
  ( 
   let red = int_of_float rgb.(0) in
   let red = if red > 255 then 255 else red in
   print_byte red;

   let green = int_of_float rgb.(1) in
   let green = if green > 255 then 255 else green in
   print_byte green;

   let blue = int_of_float rgb.(2) in
   let blue = if blue > 255 then 255 else blue in
   print_byte blue
  )
in

(*MINCAML*)let rec write_ppm_header _ =
(*NOMINCAML let write_ppm_header _ =*)
  ( 
    print_byte 80; (* 'P' *)
    print_byte (48 + 6); (* 48 = '0' *)
    print_byte 10;
    print_int size.(0);
    print_byte 32;
    print_int size.(1);
    print_byte 10;
    print_int 255;
    print_byte 10
  )
in

(**** 1茵ゃ鴻腟昭****)
(*MINCAML*)let rec scan_point scanx =
(*NOMINCAML let rec scan_point scanx =*)
  if scanx >= size.(0) then () else
  (
    (* 莎井膩激婚罔 *)
    let sscanx = (float_of_int scanx -. scan_offset.(0)) *. scan_d.(0) in
    (* 若劫: 荀襲号 *)
    vscan.(0) <- (sscanx *. cos_v.(1) +. wscan.(0));
    vscan.(1) <- (scan_sscany.(0) *. cos_v.(0) -. vp.(1));
    vscan.(2) <- (-.sscanx *. sin_v.(1) +. wscan.(2));

    (* 荀莟∞㍽≪ *)
    let metric = sqrt ((fsqr sscanx) +. scan_met1.(0)) in
    vscan.(0) <- vscan.(0) /. metric;
    vscan.(1) <- vscan.(1) /. metric;
    vscan.(2) <- vscan.(2) /. metric;

    viewpoint.(0) <- view.(0);
    viewpoint.(1) <- view.(1);
    viewpoint.(2) <- view.(2);

    (* 蚊㏍若 *)
    rgb.(0) <- 0.0;
    rgb.(1) <- 0.0;
    rgb.(2) <- 0.0;

    (* Go! *)
    raytracing 0 1.0;

    (* 贋＜阪 *)
    write_rgb ();

    (* 罨＜鴻 *)
    scan_point (scanx + 1)
  )
in

(**** 鴻ｃ ****)
(*MINCAML*)let rec scan_line scany =
(*NOMINCAML let rec scan_line scany =*)
  if scany < size.(0) then
    ( 
(*      if dbg.(0) then
        ( 
          print_string "scanning y = ";
          print_int scany;
          print_string "\n";
        )
      else ();*)
      scan_sscany.(0) <- (
        let t = (scan_offset.(0) -. 1.0 -. float_of_int scany) in
        scan_d.(0) *. t);
      (* 莎井膩㍽≪2箙*)
      scan_met1.(0) <- fsqr scan_sscany.(0) +. 40000.0;
      (* wscan 荀莎井膩с*)
      let t1 = scan_sscany.(0) *. sin_v.(0) in
      wscan.(0) <- t1 *. sin_v.(1) -. vp.(0);
      wscan.(2) <- t1 *. cos_v.(1) -. vp.(2);
      scan_point 0;
      scan_line (scany + 1)
    )
  else
    ()
in

(**** e鴻ｃ ****)
(*MINCAML*)let rec scan_start _ =
(*NOMINCAML let scan_start _ =*)
  (
    write_ppm_header ();
    let sizex = float_of_int size.(0) in
    scan_d.(0) <- 128.0 /. sizex; 
    scan_offset.(0) <- sizex /. 2.0;
    scan_line 0
  )
in

(**************** ＜潟若****************)

(*MINCAML*)let rec rt size_x size_y debug_p =
(*NOMINCAML let rt size_x size_y debug_p =*)
  ( 
    size.(0) <- size_x;
    size.(1) <- size_y;
    dbg.(0) <- debug_p;
    read_parameter ();
    scan_start ()
  )
in

rt 768 768 false
