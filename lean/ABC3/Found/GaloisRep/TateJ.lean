import ABC3.Found.GaloisRep.TateUnit

/-!
# Galois (G6) 第 99 ブロック —— **★★★★`1/j = q·(単元)`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★j 反転の形が出た

Tate 母数 `q` は、古典的には `j(q) = 1/q + 744 + ⋯` を**反転**して得る。
★形式冪級数の水準では、反転できる形かどうかが問題である。

★★`c₄ = b₂² − 24b₄ = 1 − 48a₄` の**定数項は 1**——すなわち `c₄` は**単元**である。
★★★したがって `j = c₄³/Δ` と `Δ = X·(単元)`(第 98)から

    Δ / c₄³ = X · (単元)          すなわち  1/j = q·(単元)

★★★★**これが「`q` を `j` から復元できる」形である**——
`X·(単元)` は形式的にも Hensel 的にも反転できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `constantCoeff_tateC4` | ★`c₄` の定数項は 1 |
| `tateC4_isUnit` | ★★**`c₄` は単元** |
| `tateInvJ_eq_X_mul_unit` | ★★★★**`Δ/c₄³ = X·(単元)`** |
-/

namespace ABC3.Found.GaloisRep

open PowerSeries

/-- ★`c₄ = 1 − 48 a₄`。 -/
theorem tateCurve_c4_eq : tateCurve.c₄ = 1 - 48 * tateA4 := by
  simp only [WeierstrassCurve.c₄, WeierstrassCurve.b₂, WeierstrassCurve.b₄, tateCurve]
  ring

/-- ★`c₄` の定数項は 1。 -/
theorem constantCoeff_tateC4 : PowerSeries.constantCoeff tateCurve.c₄ = 1 := by
  rw [tateCurve_c4_eq]
  simp [constantCoeff_tateA4]

/-- ★★**`c₄` は単元である**。 -/
theorem tateC4_isUnit : IsUnit tateCurve.c₄ := by
  rw [PowerSeries.isUnit_iff_constantCoeff, constantCoeff_tateC4]
  exact isUnit_one

/-- ★★★★**`Δ/c₄³ = X·(単元)`**——`1/j` が `q` の単元倍であること。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★これが Tate 母数 `q` を `j` から復元する道の形である。 -/
theorem tateInvJ_eq_X_mul_unit :
    ∃ w : PowerSeries ℤ, IsUnit w ∧ tateCurve.Δ = PowerSeries.X * w * tateCurve.c₄ ^ 3 := by
  obtain ⟨u, hu, hΔ⟩ := tateDelta_eq_X_mul_unit
  obtain ⟨c, hc⟩ := (tateC4_isUnit.pow 3)
  refine ⟨u * (↑c⁻¹ : PowerSeries ℤ), hu.mul (Units.isUnit c⁻¹), ?_⟩
  rw [hΔ, ← hc]
  rw [mul_assoc, mul_assoc]
  congr 1
  simp

/-! ## ★出典の紐付け(`.src`) -/

def tateInvJ_eq_X_mul_unit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 母数を j から復元する形)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
