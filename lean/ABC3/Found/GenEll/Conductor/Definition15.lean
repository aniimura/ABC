/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ProductFormula

/-!
# Conductor —— `[GenEll] Definition 1.5` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain
variable {F : Type*} [Field F] [NumberField F]

/-! ## ★被約化 -/

/-- **[GenEll] Definition 1.5, (ii)** の `(−)_red`、算術因子の上での形。

原文 (GenEll p.8):
> is also an effective Cartier divisor. We shall say that E is reduced if E = Ered.

★各有限素点の係数を、正なら `1`、そうでなければ `0` に潰す。
`Finsupp.mapRange` が `0 ↦ 0` を要求するので、その条件はここで満たされている。
★アルキメデス側は `0`——原文の `f_x^D` は `𝕍(F)^non` に台を持つからである。 -/
noncomputable def ADivRed (a : ADiv F) : ADiv F :=
  (Finsupp.mapRange (fun n : ℤ => if 0 < n then (1 : ℤ) else 0) (by simp) a.fin, 0)

/-! ## ★出典の紐付け(`.src`) -/

/-- ★**条を明示する**(2026-08-17 夕に修正)。

`ADivRed` は `ADiv F`(= `Spec 𝓞_F` 上の算術因子)の被約化であり、
原文 (ii) が言う**一般の正規ネーター scheme `Z`** の場合ではない。
★`Spec 𝓞_F` は Dedekind なので非零イデアルはすべて可逆であり、
**正則性も Auslander–Buchsbaum も要らない**——だから先に取れた。

★★一般の場合は `RadicalCartier.lean` にあり、そこには
「各茎が UFD」という仮定が残っている(正則 ⟹ UFD が mathlib に無いため)。
★**ラベルが条を書いていないと、一般の場合まで取れたように読める。** -/
def ADivRed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (ii)(Spec 𝓞_F の場合のみ——一般の Z は RadicalCartier.lean)",
    sectionId := "genell-def-1-5" }

end ABC3.Found.GenEll
