open MiniMLRuntime;;

(**************** 違若違絎ｈ ****************)

(* ⒢吾若鐚紊0*)
let objects = 
  let dummy = Array.make 0 0.0 in
  Array.make (60) 
    (0, 0, 0, 0, 
     dummy, dummy,
     false, dummy, dummy,
     dummy)

(* [| x荵吾莎井膩, y荵吾莎井膩|] *)
let size = Array.make 2 128

(* 絎с: 阪*)
let dbg = Array.make 1 true
(* Screen 婚罔*)
let screen = Array.make 3 0.0
(* 荀婚罔(offset ⒢) *)
let vp = Array.make 3 0.0
(* 荀婚罔(screen 篏舟 offset ) *)
let view = Array.make 3 0.0
(* 劫 () *)
let light = Array.make 3 0.0
(* 鴻⒢潟拶劫: 筝∽$т*)
let cos_v = Array.make 2 0.0
let sin_v = Array.make 2 0.0
(* ♂ゃゃ綣桁墾 (罔=255) *)
let beam = Array.make 1 255.0
(* AND 若*)
let and_net = Array.make 50 (Array.make 1 (-1))
(* OR 若*)
let or_net = Array.make 1 (Array.make 1 (and_net.(0)))

(* reader *)
let temp = Array.make 14 0.0 (* read_nth_object 篏キ紊 *)
let cs_temp = Array.make 16 0.0

(* solver *)
(**** Callee 拭㏍紊 ****)
(* 篋ょ t $ *)
let solver_dist = Array.make 1 0.0

(* 鴻ｃ*)
let vscan = Array.make 3 0.0
(* 篋ょ㍾剛茵с劫 *)
let intsec_rectside = Array.make 1 0
(* 肴篋ょ絨t *)
let tmin = Array.make 1 (1000000000.0)
(* 篋ょ婚罔*)
let crashed_point = Array.make 3 0.0
(* 茵⒢吾 *)
let crashed_object = Array.make 1 0
(* 1ゃ AND 若㍾篋❼ *)
let end_flag = Array.make 1 false
(* 若*)
let viewpoint = Array.make 3 0.0
(* 羈 *)
let nvector = Array.make 3 0.0
(* 鴻⒢割㍾ *)
let rgb = Array.make 3 0.0
(* 篋ょ㍽ *)
let texture_color = Array.make 3 0.0

(* ⒢吾筝㊤ 鴻荀 *)
let solver_w_vec = Array.make 3 0.0

(* check_all_inside 違*)
let chkinside_p = Array.make 3 0.0

(* is_outside (筝㊤綏)  *)
let isoutside_q = Array.make 3 0.0

(* 違若冴若*)
(* nvector *)
let nvector_w = Array.make 3 0.0

(* main *)
let scan_d = Array.make 1 0.0
let scan_offset = Array.make 1 0.0
let scan_sscany = Array.make 1 0.0
let scan_met1 = Array.make 1 0.0
let wscan = Array.make 3 0.0
