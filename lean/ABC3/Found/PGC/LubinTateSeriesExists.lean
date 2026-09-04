import ABC3.Found.PGC.LubinTateFormalGroupLaw

/-!
# Lubin-Tate 級数 `f` は存在する(`sorry` 無し)

`Found/PGC/LubinTateFormalGroupLaw.lean`(および、それに続く一連の
`Found/PGC/LubinTate*.lean`)は、`f≡πX(\mathrm{mod}\deg2)`・
`f≡X^q(\mathrm{mod}\ π)` を満たす冪級数 `f` を**仮説として**取り、
その `f` から形式群法則 `F_f`・`𝒪_K` 作用 `[a]_f` 等を構成してきた。
本ファイルは、そのような `f` が実際に**存在する**ことを示す——
古典的には `f(X):=πX+X^q` という最も単純な多項式で足りる。

これにより、これまでの `Found/PGC/LubinTateAction*.lean` 一連の定理群
(`F_f` の存在・可換律・結合律・𝒪_K 作用への拡張とその環準同型性・
Weierstrass 標準分解・`ψ_n` の既約性等)が**空虚な仮説の上の議論ではない**
ことが確定する——`f` の存在は付値環の基本的な性質(極大イデアルが `π`
で生成されること)だけから直ちに従う、ごく初等的な事実。
-/

namespace ABC3.Found.PGC

/-- ★★★★★★★★★**Lubin-Tate 級数 `f` は存在する**——`f(X):=πX+X^q`
(`q:=pp^ff`)が `f≡πX(\mathrm{mod}\deg2)`・`f≡X^q(\mathrm{mod}\ π)`
を満たす。`X^q` の次数0・次数1の係数が0であること(`q≥2`、剰余体が
真の体で位数`≥2`であることから)と、`π` が剰余への還元で0になること
(`π∈maximalIdeal`)だけから出る——新しい構成は一切不要だった。 -/
theorem exists_lubinTateSeries {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) :
    ∃ f : PowerSeries A, PowerSeries.coeff 0 f = 0 ∧ PowerSeries.coeff 1 f = π ∧
      PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff) := by
  have hcard2 : 1 < pp ^ ff := hq ▸ Fintype.one_lt_card
  refine ⟨PowerSeries.C π * PowerSeries.X + PowerSeries.X ^ (pp ^ ff), ?_, ?_, ?_⟩
  · show PowerSeries.coeff 0 (PowerSeries.C π * PowerSeries.X) +
        PowerSeries.coeff 0 (PowerSeries.X ^ (pp ^ ff) : PowerSeries A) = 0
    have h1 : PowerSeries.coeff 0 (PowerSeries.C π * PowerSeries.X : PowerSeries A) = 0 := by simp
    have h2 : PowerSeries.coeff 0 (PowerSeries.X ^ (pp ^ ff) : PowerSeries A) = 0 := by
      rw [PowerSeries.coeff_zero_eq_constantCoeff, map_pow, PowerSeries.constantCoeff_X]
      exact zero_pow (by omega)
    rw [h1, h2, add_zero]
  · show PowerSeries.coeff 1 (PowerSeries.C π * PowerSeries.X) +
        PowerSeries.coeff 1 (PowerSeries.X ^ (pp ^ ff) : PowerSeries A) = π
    have h1 : PowerSeries.coeff 1 (PowerSeries.C π * PowerSeries.X : PowerSeries A) = π := by simp
    have h2 : PowerSeries.coeff 1 (PowerSeries.X ^ (pp ^ ff) : PowerSeries A) = 0 := by
      rw [PowerSeries.coeff_X_pow]
      simp
      omega
    rw [h1, h2, add_zero]
  · have hres : IsLocalRing.residue A π = 0 := by
      rw [IsLocalRing.residue_eq_zero_iff, hπmax]
      exact Ideal.subset_span rfl
    rw [map_add, map_mul, map_pow]
    simp [hres]

end ABC3.Found.PGC
