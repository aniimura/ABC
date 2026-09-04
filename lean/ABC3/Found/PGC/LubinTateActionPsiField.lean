import ABC3.Found.PGC.LubinTateActionPsi

/-!
# `ψ_n` の根を添加した体は次数 `q^n-q^{n-1}`(`sorry` 無し)

`Found/PGC/LubinTateActionPsi.lean::irreducible_iteratedLubinTatePsi`
(`ψ_n` は `A` 上既約)を、フラクション体 `FractionRing A` 上での既約性
(Gauss の補題)へ橋渡しし、`ψ_n` の根を1つ添加した体の拡大次数が
古典的な値 `q^n-q^{n-1}` に一致することを示す。

## 逸脱の記録(CLAUDE.md 逸脱)

Gauss の補題(`Polynomial.IsPrimitive.irreducible_iff_irreducible_
map_fraction_map`)は `[Nonempty (NormalizedGCDMonoid R)]`(または
同値な `UniqueFactorizationMonoid R`)を要求する。古典的な Lubin-Tate
理論では `A=𝒪_K` は完備離散付値環(DVR)であり、これは自動的に
成り立つ性質だが、この節目のこれまでの議論(`[CommRing A][IsLocalRing
A][IsDomain A][IsAdicComplete (maximalIdeal A) A]`)には**含まれて
いなかった**。ここで**新たに `[UniqueFactorizationMonoid A]` を追加
する**——既存の定理(`irreducible_iteratedLubinTatePsi` 等)には一切
触れず、この文書内の新しい定理だけに使うので、後続の証明への影響は
無い。場所: 本ファイル全体。理由: 古典的設定(`A` は DVR)で自動的に
成り立つ性質を明示化しただけであり、原典の意図(𝒪_K は DVR)に忠実。
-/

namespace ABC3.Found.PGC

/-- ★★★★★★★★**`ψ_n` は `FractionRing A` 上でも既約**——Gauss の補題
(`Polynomial.IsPrimitive.irreducible_iff_irreducible_map_fraction_map`、
`ψ_n` がモニックなので `IsPrimitive`)で `A` 上の既約性
(`irreducible_iteratedLubinTatePsi`)を橋渡しする。 -/
theorem irreducible_iteratedLubinTatePsi_map_fractionRing {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A] [UniqueFactorizationMonoid A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    Irreducible (Polynomial.map (algebraMap A (FractionRing A))
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)) := by
  have hmonic := (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic
  exact (Polynomial.IsPrimitive.irreducible_iff_irreducible_map_fraction_map hmonic.isPrimitive).mp
    (irreducible_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)

/-- ★★★★★★★★★**`ψ_n` の根を1つ添加した体の `FractionRing A` 上の
拡大次数は `q^n-q^{n-1}`**——古典的な Lubin-Tate 理論の中心的な結論
(原始 `π^n`-捩れ点を1つ添加した拡大の次数)を体の言葉で述べたもの。
`AdjoinRoot` の `PowerBasis`(既約多項式の根を添加した体の標準基底)の
次元が多項式の次数に一致するという一般事実
(`AdjoinRoot.powerBasis_dim`)と、`ψ_n` の次数公式
(`natDegree_iteratedLubinTatePsi`)を組み合わせるだけ。 -/
theorem finrank_adjoinRoot_iteratedLubinTatePsi {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A] [UniqueFactorizationMonoid A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    Module.finrank (FractionRing A) (AdjoinRoot (Polynomial.map (algebraMap A (FractionRing A))
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn))) =
      (pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) := by
  have hne0 : Polynomial.map (algebraMap A (FractionRing A))
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn) ≠ 0 :=
    (irreducible_iteratedLubinTatePsi_map_fractionRing hq hπmax hπne0 f hf0 hf1 hf n hn).ne_zero
  rw [(AdjoinRoot.powerBasis hne0).finrank, AdjoinRoot.powerBasis_dim hne0,
    (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic.natDegree_map,
    natDegree_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn]

end ABC3.Found.PGC
