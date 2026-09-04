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

/-! ### 分離性: `ψ_n` の根は代数閉体で相異なる `q^n-q^{n-1}` 個(混標数の場合)

`[CharZero A]` を追加する(古典的な Lubin-Tate 理論の**混標数**の場合
——`K` は `ℚ_p` の有限次拡大——に対応する仮定、CLAUDE.md 逸脱の記録:
既存の定理には一切触れず、この節の新しい定理だけに使う)。標数0の体
では既約多項式は自動的に分離的(`Irreducible.separable`)なので、
`ψ_n` の根が**相異なる**(重複度が無い)ことまで結論できる——前回
`card_iteratedLubinTateDistinguishedRoots` 等で「重複度込み」としか
言えなかったギャップが、混標数の場合には解消される。 -/

/-- `ψ_n` は `FractionRing A` 上で分離的——標数0では既約多項式は
自動的に分離的(`Irreducible.separable`)。 -/
theorem separable_iteratedLubinTatePsi_map_fractionRing {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A] [UniqueFactorizationMonoid A] [CharZero A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    (Polynomial.map (algebraMap A (FractionRing A))
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)).Separable :=
  (irreducible_iteratedLubinTatePsi_map_fractionRing hq hπmax hπne0 f hf0 hf1 hf n hn).separable

/-- `ψ_n` の根(`iteratedLubinTateAlgClosure A` の中、重複度込み)。 -/
noncomputable def iteratedLubinTatePsiRoots {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    Multiset (iteratedLubinTateAlgClosure A) :=
  (Polynomial.map (algebraMap A (iteratedLubinTateAlgClosure A))
    (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)).roots

/-- `ψ_n` の根は重複度込みでちょうど `q^n-q^{n-1}` 個——
`card_iteratedLubinTateDistinguishedRoots` と全く同じ議論。 -/
theorem card_iteratedLubinTatePsiRoots {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    Multiset.card (iteratedLubinTatePsiRoots hq hπmax hπne0 f hf0 hf1 hf n hn) =
      (pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) := by
  rw [iteratedLubinTatePsiRoots, IsAlgClosed.card_roots_eq_natDegree,
    Polynomial.Monic.natDegree_map (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic
      (algebraMap A (iteratedLubinTateAlgClosure A)),
    natDegree_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn]

/-- ★★★★★★★★★★**`ψ_n` の根は互いに相異なる**(混標数、`[CharZero A]`
のとき)——`ψ_n` の `FractionRing A` 上の分離性
(`separable_iteratedLubinTatePsi_map_fractionRing`)を
`iteratedLubinTateAlgClosure A` へ持ち上げ(`Polynomial.Separable.map`、
代入の合成 `IsScalarTower.algebraMap_eq`)、`Polynomial.nodup_roots`
を適用する。`card_iteratedLubinTatePsiRoots`(個数は`q^n-q^{n-1}`)と
合わせて、混標数の場合には**原始 `π^n`-捩れ点が真に `q^n-q^{n-1}` 個の
相異なる元からなる**という古典的な Lubin-Tate 理論の帰結が完成する。 -/
theorem nodup_iteratedLubinTatePsiRoots {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A] [UniqueFactorizationMonoid A] [CharZero A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (n : ℕ) (hn : 1 ≤ n) :
    (iteratedLubinTatePsiRoots hq hπmax hπne0 f hf0 hf1 hf n hn).Nodup := by
  show (Polynomial.map (algebraMap A (iteratedLubinTateAlgClosure A))
    (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)).roots.Nodup
  apply Polynomial.nodup_roots
  rw [IsScalarTower.algebraMap_eq A (FractionRing A) (iteratedLubinTateAlgClosure A),
    ← Polynomial.map_map]
  exact (separable_iteratedLubinTatePsi_map_fractionRing hq hπmax hπne0 f hf0 hf1 hf n hn).map

end ABC3.Found.PGC
