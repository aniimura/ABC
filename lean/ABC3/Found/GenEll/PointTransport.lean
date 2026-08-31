/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateVarChange
import ABC3.Found.GenEll.PointVariableChange

/-!
# 第 966 ブロック —— **★★★★★★★★★★★★捉れ点を運ぶ**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か

`minDeltaExp_eq_mul_of_torsion`（第 965）は Tate モデルの上の点 `P` で

    `l • P = 0`  かつ  `P ≠ 0`

を受ける。☆一方 `Lemma 3.5` が持っているのは **`L` の上の点 `Q`**（`addOrderOf Q = l`）である。

★その隔たりは 2 つの写像でできている:

* `rhPoint`——体準同型 `L → Lv` による輸送
* `vcPoint`——変数変換 `C` による輸送（`E ⊗ Lv` から Tate モデルへ）

☆どちらも「`l • Q = 0` を保つ」ことと「`0` でないものを `0` にしない」ことだけが要る。
★**位数がちょうど保たれることまでは要らない**——それが証明を軽くする。

## ★配管の注意（実測）

`rhPoint_nsmul` などの `•` は `open scoped Classical` の下で
`Point` の加法群を取っている。☆こちらが `[DecidableEq F]` を束縛すると
**別の `SMul` インスタンス**になって `rw` が当たらない。
★`open scoped Classical in` を付けて `DecidableEq` は束縛しないこと。

| 定理 | 内容 |
|---|---|
| `nsmul_eq_zero_of_addOrderOf` | ★`addOrderOf Q = l` なら `l • Q = 0` |
| `ne_zero_of_addOrderOf_prime` | ★`l` が素なら `Q ≠ 0` |
| `rhPoint_nsmul_eq_zero` / `rhPoint_ne_zero` | ★★★★体拡大で運ぶ |
| `vcPoint_nsmul_eq_zero` / `vcPoint_ne_zero` | ★★★★変数変換で運ぶ |
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep

/-! ## ★位数から 2 つの事実へ -/

/-- ★**`addOrderOf Q = l` なら `l • Q = 0`**。 -/
theorem nsmul_eq_zero_of_addOrderOf {G : Type} [AddGroup G] {Q : G} {l : ℕ}
    (hQ : addOrderOf Q = l) : l • Q = 0 := by
  rw [← hQ]; exact addOrderOf_nsmul_eq_zero Q

/-- ★**位数が素数なら `Q ≠ 0`**。 -/
theorem ne_zero_of_addOrderOf_prime {G : Type} [AddGroup G] {Q : G} {l : ℕ} (hl : l.Prime)
    (hQ : addOrderOf Q = l) : Q ≠ 0 := by
  intro hc
  rw [hc, addOrderOf_zero] at hQ
  exact hl.one_lt.ne hQ

/-! ## ★★★★体拡大で運ぶ -/

open scoped Classical in
/-- ★★★★**体拡大で `l`-捉れは `l`-捉れのまま**。 -/
theorem rhPoint_nsmul_eq_zero {F K : Type} [Field F] [Field K]
    (f : F →+* K) (W : WeierstrassCurve F) [W.IsElliptic] [(W.map f).IsElliptic]
    {Q : W.toAffine.Point} {n : ℕ} (h : n • Q = 0) : n • rhPoint f W Q = 0 := by
  rw [← rhPoint_nsmul, h, rhPoint_zero]

/-- ★★★★**体拡大で `O` でない点は `O` にならない**——`rhPoint` は単射だから。 -/
theorem rhPoint_ne_zero {F K : Type} [Field F] [Field K] (f : F →+* K)
    (W : WeierstrassCurve F) {Q : W.toAffine.Point} (hQ : Q ≠ 0) :
    rhPoint f W Q ≠ 0 := by
  intro hc
  exact hQ (rhPoint_injective f W (by rw [hc, rhPoint_zero]))

/-! ## ★★★★変数変換で運ぶ -/

open scoped Classical in
/-- ★★★★**変数変換で `l`-捉れは `l`-捉れのまま**。 -/
theorem vcPoint_nsmul_eq_zero {F : Type} [Field F]
    (C : WeierstrassCurve.VariableChange F) (W : WeierstrassCurve F)
    [W.IsElliptic] [(C • W).IsElliptic]
    {Q : W.toAffine.Point} {n : ℕ} (h : n • Q = 0) : n • vcPoint C W Q = 0 := by
  rw [← vcPoint_nsmul, h, ABC3.Found.GenEll.vcPoint_zero]

open scoped Classical in
/-- ★★★★**変数変換で `O` でない点は `O` にならない**。 -/
theorem vcPoint_ne_zero {F : Type} [Field F]
    (C : WeierstrassCurve.VariableChange F)
    (W : WeierstrassCurve F) {Q : W.toAffine.Point} (hQ : Q ≠ 0) :
    vcPoint C W Q ≠ 0 := by
  intro hc
  refine hQ ((vcPointAddEquiv W C).injective ?_)
  rw [show (vcPointAddEquiv W C) Q = vcPoint C W Q from rfl, hc, map_zero]

/-! ## ★出典の紐付け(`.src`) -/

def rhPoint_nsmul_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(体拡大で l-捉れは l-捉れのまま。★無条件)",
    sectionId := "genell-lemma-3-5" }

def vcPoint_nsmul_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(変数変換で l-捉れは l-捉れのまま。★無条件)",
    sectionId := "genell-lemma-3-5" }

def rhPoint_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(体拡大で O でない点は O にならない。★無条件)",
    sectionId := "genell-lemma-3-5" }

def vcPoint_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(変数変換で O でない点は O にならない。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
