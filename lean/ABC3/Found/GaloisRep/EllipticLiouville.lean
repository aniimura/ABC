import ABC3.Found.GaloisRep.TateEquation

/-!
# Galois (G6) 第 241 ブロック —— **★★★★★★楕円関数の Liouville**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★葉 (c) の道を測った

葉 (b) が閉じた(第 240)。次は**葉 (c)——写像 `Kˣ/q^ℤ → E(K)` が準同型であること**である。
古典的にはこれは**三点の共線性**である:

    u₁ u₂ u₃ = 1  ⟹  (X(u₁),Y(u₁)), (X(u₂),Y(u₂)), (X(u₃),Y(u₃)) は一直線上

★★★在庫を引いた結果:**mathlib に `℘` の加法定理は無い**。
`Analysis/SpecialFunctions/Elliptic/Weierstrass.lean` は `derivWeierstrassP_sq`
(`℘'² = 4℘³ − g₂℘ − g₃`)で終わっており、加法定理も「楕円関数は定数」も無い。

★★★★したがって道は次の順になる:

| 段 | 内容 | 状態 |
|---|---|---|
| 1 | 楕円関数の Liouville(二重周期な整関数は定数) | ★**本ブロック** |
| 2 | `℘` の加法定理(共線性の形) | 未 |
| 3 | 解析側の共線性を第 221-240 の機構で形式側へ移す | 未 |

★段 3 は**葉 (b) で作った機構がそのまま使える**——変数が `(u₁, u₂, q)` の 3 つに増えるだけで、
万有な環・特殊化・係数の消滅の流れは同じである。

## ★★★★★★段 1 —— 楕円関数の Liouville

道具はすべて mathlib に在った:

| 部品 | mathlib |
|---|---|
| Liouville | `Differentiable.apply_eq_apply_of_bounded` |
| 基本領域の存在 | `ZSpan.exist_unique_vadd_mem_fundamentalDomain` |
| 基本領域の有界性 | `ZSpan.fundamentalDomain_isBounded` |
| 有界閉集合はコンパクト | `Bornology.IsBounded.isCompact_closure`(ℂ は proper) |

★**値域は基本領域の閉包の像に含まれる**——周期性で任意の点を基本領域に移せるから。
そこは連続像なのでコンパクト、したがって有界。あとは Liouville。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `eq_of_periodic_differentiable` | ★★★★★★**楕円関数の Liouville** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real PeriodPair

/-- ★★★★★★**楕円関数の Liouville**——格子で二重周期的な整関数は定数である。

★値域は基本領域の閉包の像に含まれる(周期性)。そこは連続像なのでコンパクト、
したがって有界。あとは `Differentiable.apply_eq_apply_of_bounded`。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem eq_of_periodic_differentiable (L : PeriodPair) (f : ℂ → ℂ)
    (hdiff : Differentiable ℂ f)
    (hper : ∀ l ∈ L.lattice, ∀ z : ℂ, f (l + z) = f z) :
    ∀ z w : ℂ, f z = f w := by
  have hsub : Set.range f ⊆ f '' (closure (ZSpan.fundamentalDomain L.basis)) := by
    rintro y ⟨x, rfl⟩
    obtain ⟨v, hv, -⟩ := ZSpan.exist_unique_vadd_mem_fundamentalDomain L.basis x
    have hvmem : (v : ℂ) ∈ L.lattice := by
      rw [L.lattice_eq_span_range_basis]
      exact v.2
    exact ⟨(v : ℂ) + x, subset_closure hv, hper (v : ℂ) hvmem x⟩
  have hbd : Bornology.IsBounded (Set.range f) := by
    refine Bornology.IsBounded.subset ?_ hsub
    refine (IsCompact.image ?_ hdiff.continuous).isBounded
    exact (ZSpan.fundamentalDomain_isBounded L.basis).isCompact_closure
  exact fun z w => hdiff.apply_eq_apply_of_bounded hbd z w

/-! ## ★出典の紐付け(`.src`) -/

def eq_of_periodic_differentiable.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——楕円関数の Liouville)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
