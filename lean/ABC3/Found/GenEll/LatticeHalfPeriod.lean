import ABC3.Meta.Claim
import Mathlib.Analysis.SpecialFunctions.Elliptic.Weierstrass

/-!
# GenEll 第 332 ブロック —— **★★★★★半周期と 3 次式の根**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★これは G8 の一意化の連鎖の最初の段である

`Skeleton/GenEll/Uniformization.lean` の `exists_periodPair`(任意の `E/ℂ` が `ℂ/Λ`)は
G8 の `htFalt` を固定するための最下流の葉である。★その手前に

    (i)  束の判別式 `g₂³ − 27g₃² ≠ 0`
    (ii) `ℂ/Λ ≅ E(ℂ)` の群同型
    (iii) `j` の全射性

の 3 段があり、**mathlib はいずれも持たない**(2026-08-26 実測)。
★★本ブロックは (i) の最初の段——**半周期での `℘'` の消失と 3 次式の根**——を取る。

## ★★★★★中身

`℘'` は**奇関数**(`derivWeierstrassP_neg`)かつ**周期的**(`derivWeierstrassP_add_coe`)なので、
`w ∈ Λ` に対し

    ℘'(w/2) = ℘'(-w/2 + w) = ℘'(-w/2) = -℘'(w/2)

より `℘'(w/2) = 0`。★★★これを `℘'² = 4℘³ − g₂℘ − g₃`(`derivWeierstrassP_sq`)に入れると
`℘(w/2)` は 3 次式 `4x³ − g₂x − g₃` の**根**である。

★★★★半周期は 3 つある——`ω₁/2`・`ω₂/2`・`(ω₁+ω₂)/2`。
前 2 つが束に入らないことは mathlib にあり(`ω₁_div_two_notMem_lattice` 等)、
3 つ目は `mul_ω₁_add_mul_ω₂_mem_lattice`(有理係数の判定)から出る。

## ★★残り((i) の続き)

★3 つの根が**相異なる**ことが要る。それには「楕円関数の零点と極の個数は等しい」
——基本平行四辺形の周での偏角の原理——が要り、mathlib には無い。
★★`℘ − e` は 0 で 2 位の極を持つので零点はちょうど 2 個、
一方 `e₁ = e₂` なら `ω₁/2` と `ω₂/2` で各 2 位、計 4 個になって矛盾する、という筋である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `derivWeierstrassP_half` | ★★★★**半周期で `℘' = 0`** |
| `isRoot_cubic_weierstrassP_half` | ★★★★★**`℘(w/2)` は `4x³ − g₂x − g₃` の根** |
| `add_div_two_notMem_lattice` | ★★★第 3 の半周期も束に入らない |
| `isRoot_cubic_e₁`・`_e₂`・`_e₃` | ★★★★3 つの根 |
-/

namespace ABC3.Found.GenEll

open PeriodPair

/-! ## ★★★★半周期で `℘'` が消える -/

/-- ★★★★**半周期では `℘'` が消える**——奇関数性と周期性から。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any -/
theorem derivWeierstrassP_half (L : PeriodPair) (w : ℂ) (hw : w ∈ L.lattice) :
    L.derivWeierstrassP (w / 2) = 0 := by
  have h1 : L.derivWeierstrassP (-(w / 2) + w) = L.derivWeierstrassP (-(w / 2)) :=
    L.derivWeierstrassP_add_coe (-(w / 2)) ⟨w, hw⟩
  rw [show -(w / 2) + w = w / 2 by ring, L.derivWeierstrassP_neg] at h1
  linear_combination h1 / 2

/-- ★★★★★**半周期での `℘` の値は 3 次式 `4x³ − g₂x − g₃` の根である**。 -/
theorem isRoot_cubic_weierstrassP_half (L : PeriodPair) (w : ℂ) (hw : w ∈ L.lattice)
    (hw2 : w / 2 ∉ L.lattice) :
    4 * L.weierstrassP (w / 2) ^ 3 - L.g₂ * L.weierstrassP (w / 2) - L.g₃ = 0 := by
  have h := L.derivWeierstrassP_sq (w / 2) hw2
  rw [derivWeierstrassP_half L w hw] at h
  linear_combination -h

/-! ## ★★★3 つの半周期 -/

/-- ★★★第 3 の半周期 `(ω₁+ω₂)/2` も束に入らない。 -/
theorem add_div_two_notMem_lattice (L : PeriodPair) : (L.ω₁ + L.ω₂) / 2 ∉ L.lattice := by
  intro h
  have hq : ((1 / 2 : ℚ) : ℂ) * L.ω₁ + ((1 / 2 : ℚ) : ℂ) * L.ω₂ ∈ L.lattice := by
    convert h using 1
    push_cast
    ring
  rw [PeriodPair.mul_ω₁_add_mul_ω₂_mem_lattice] at hq
  norm_num at hq

/-- ★★★★`e₁ := ℘(ω₁/2)` は 3 次式の根。 -/
theorem isRoot_cubic_e₁ (L : PeriodPair) :
    4 * L.weierstrassP (L.ω₁ / 2) ^ 3 - L.g₂ * L.weierstrassP (L.ω₁ / 2) - L.g₃ = 0 :=
  isRoot_cubic_weierstrassP_half L L.ω₁ L.ω₁_mem_lattice L.ω₁_div_two_notMem_lattice

/-- ★★★★`e₂ := ℘(ω₂/2)` は 3 次式の根。 -/
theorem isRoot_cubic_e₂ (L : PeriodPair) :
    4 * L.weierstrassP (L.ω₂ / 2) ^ 3 - L.g₂ * L.weierstrassP (L.ω₂ / 2) - L.g₃ = 0 :=
  isRoot_cubic_weierstrassP_half L L.ω₂ L.ω₂_mem_lattice L.ω₂_div_two_notMem_lattice

/-- ★★★★`e₃ := ℘((ω₁+ω₂)/2)` は 3 次式の根。 -/
theorem isRoot_cubic_e₃ (L : PeriodPair) :
    4 * L.weierstrassP ((L.ω₁ + L.ω₂) / 2) ^ 3
      - L.g₂ * L.weierstrassP ((L.ω₁ + L.ω₂) / 2) - L.g₃ = 0 :=
  isRoot_cubic_weierstrassP_half L (L.ω₁ + L.ω₂)
    (L.lattice.add_mem L.ω₁_mem_lattice L.ω₂_mem_lattice) (add_div_two_notMem_lattice L)

/-! ## ★出典の紐付け(`.src`) -/

def derivWeierstrassP_half.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def isRoot_cubic_weierstrassP_half.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
