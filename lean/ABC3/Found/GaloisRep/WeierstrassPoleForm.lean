import ABC3.Found.GaloisRep.WeierstrassPole

/-!
# Galois (G6) 第 244 ブロック —— **★★★★★極が相殺する鍵**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★相殺は「位数を数える」のではなく「因数分解」で出る

加法定理の証明で使う関数

    F(z) = ℘(z+w) + ℘(z) + ℘(w) − (1/4)·((℘'(z)−℘'(w))/(℘(z)−℘(w)))²

が格子点 `l₀` のまわりで正則であることを示したい。`t := z − l₀` として第 243 より

    ℘(z) = 1/t² + E(z),   ℘'(z) = −2/t³ + D(z)      (`E, D` は `l₀` で解析的)

とおくと、`c := ℘(w)`、`c' := ℘'(w)`、`α := D − c'`、`β := E − c` として

    (℘'(z)−c')/(℘(z)−c) = (1/t)·N/M,   N = −2 + t³α,   M = 1 + t²β

となる。ここで

    (1/4)·(N/M)²·(1/t²) − 1/t² = (N² − 4M²)/(4t²M²)

であり、★★**`N² − 4M² = t²·R`(`R` は解析的)と因数分解できる**(`pole_cancel_factor`)。
したがって上の式は `R/(4M²)` となり **`l₀` で解析的**である(`M(l₀) = 1 ≠ 0`)。

★★★**位数を数える必要はない**——`ring` で確かめられる恒等式ひとつで済む。
`AnalyticAt.order` を持ち出さずに相殺が出るのがこの形の利点である。

## ★★見積もりの再訂正

§9-556 で 15–35 ブロックとしたが、この因数分解で相殺の議論が機械化できるので
**8–20 ブロック**が見込みである。悪い点は 3 種類(`z ∈ L`、`z ∈ L−w`、`z ≡ w`)で、
どれも同じ形の因数分解で処理できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_pole_form` | ★★★★格子点のまわりの `℘`・`℘'` の形 |
| `differentiableOn_weierstrassP_compl` 他 | ★格子の外での微分可能性(向きを揃えた形) |
| `weierstrassP_lattice_add` 他 | ★★周期性(`l + z` の向き) |
| `pole_cancel_factor` | ★★★★★**極が相殺する鍵の因数分解** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real PeriodPair

/-! ## ★★★★格子点のまわりの形 -/

/-- ★★★★**格子点のまわりの `℘`・`℘'` の形**——主要部 + `l₀` で解析的な残り。 -/
theorem exists_pole_form (L : PeriodPair) (l₀ : L.lattice) :
    ∃ E D : ℂ → ℂ, AnalyticAt ℂ E (l₀ : ℂ) ∧ AnalyticAt ℂ D (l₀ : ℂ)
      ∧ (∀ z : ℂ, L.weierstrassP z = 1 / (z - (l₀ : ℂ)) ^ 2 + E z)
      ∧ (∀ z : ℂ, L.derivWeierstrassP z = -2 / (z - (l₀ : ℂ)) ^ 3 + D z) := by
  refine ⟨fun z => L.weierstrassPExcept (l₀ : ℂ) z - 1 / (l₀ : ℂ) ^ 2,
    fun z => L.derivWeierstrassPExcept (l₀ : ℂ) z, ?_, ?_, ?_, ?_⟩
  · exact (L.analyticAt_weierstrassPExcept (l₀ : ℂ)).sub analyticAt_const
  · exact L.analyticAt_derivWeierstrassPExcept (l₀ : ℂ)
  · intro z
    exact weierstrassP_principal L l₀ z
  · intro z
    exact derivWeierstrassP_principal L l₀ z

/-! ## ★微分可能性と周期性(向きを揃えた形) -/

theorem differentiableOn_weierstrassP_compl (L : PeriodPair) :
    DifferentiableOn ℂ L.weierstrassP ((L.lattice : Set ℂ)ᶜ) :=
  L.differentiableOn_weierstrassP

theorem differentiableOn_derivWeierstrassP_compl (L : PeriodPair) :
    DifferentiableOn ℂ L.derivWeierstrassP ((L.lattice : Set ℂ)ᶜ) :=
  L.differentiableOn_derivWeierstrassP

theorem weierstrassP_lattice_add (L : PeriodPair) (l : ℂ) (hl : l ∈ L.lattice) (z : ℂ) :
    L.weierstrassP (l + z) = L.weierstrassP z := by
  rw [show l + z = z + l by ring]
  exact L.weierstrassP_add_coe z ⟨l, hl⟩

theorem derivWeierstrassP_lattice_add (L : PeriodPair) (l : ℂ) (hl : l ∈ L.lattice) (z : ℂ) :
    L.derivWeierstrassP (l + z) = L.derivWeierstrassP z := by
  rw [show l + z = z + l by ring]
  exact L.derivWeierstrassP_add_coe z ⟨l, hl⟩

/-! ## ★★★★★極が相殺する鍵 -/

/-- ★★★★★**極が相殺する鍵の因数分解**——`N² − 4M²` は `t²` で割り切れる。

`N = −2 + t³α`(`℘'` の側)、`M = 1 + t²β`(`℘` の側)。
★**位数を数える必要はない**——`ring` で確かめられる恒等式ひとつで済む。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem pole_cancel_factor (t α β : ℂ) :
    (-2 + t ^ 3 * α) ^ 2 - 4 * (1 + t ^ 2 * β) ^ 2
      = t ^ 2 * (-4 * t * α + t ^ 4 * α ^ 2 - 8 * β - 4 * t ^ 2 * β ^ 2) := by
  ring

/-! ## ★出典の紐付け(`.src`) -/

def exists_pole_form.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——格子点のまわりの ℘・℘' の形)",
    sectionId := "genell-def-3-3" }

def pole_cancel_factor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——極が相殺する鍵の因数分解)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
