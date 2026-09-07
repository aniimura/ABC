import ABC3.Found.PGC.InertiaReduction
import ABC3.Found.PGC.FixedRingBaseAlgebra
import ABC3.Found.PGC.SubgroupCorrespondenceConstruction

/-!
# 段データの組み立て —— 各有限段の `Gal(L/L₀)^v` を作り、`Γ_K` へ引き戻す

`Found/PGC/AbsGalRamificationFiltration.lean`(Y19)は
`Interface/PGC/LocalFieldData.lean` の `RamificationFiltration p` を
「各 `K` ごとの段データ `StageFiltration K.absGal`」に還元した。
`Found/PGC/InertiaReduction.lean`(Y19b+c)は、一般の有限 Galois `L/K` を
不分岐部分 `L₀ = L ⊓ K^ur` へ底を移すことで Yoshida §6.1 の設定に落とせること
(`Gal(L/L₀) ≤ G_0`、したがって `G_0(L/L₀) = ⊤`)を示した。

**本ファイルはその段データの中身、すなわち各有限段の `Gal(L/L₀)^v` を実際に作る。**

## ★★何が埋まったか / 何が埋まっていないか(正直な線引き)

**埋まったもの**

* §1 **抽象核**(純可換環 + 群作用。分岐・付値・Galois の語彙が 1 つも出てこない):
  - `exists_fixed_sub_mem_maximalIdeal` / `exists_sub_mem_maximalIdeal_fixedRing_of_le_inertia`
    —— ★★**Teichmüller 代表元の一意性から「`H` が剰余体に自明に作用するなら
    固定環 `B^H` の剰余体は `B` の剰余体に等しい」**。
    これが原典 §6.1 の "totally ramified"(`hres` / `hresA`)の中身である。
  - `smul_mem_pow_maximalIdeal` / `normal_map_subtype_lowerRamificationGroup`
    —— `H ⊴ G` なら `G_n(H)` の `G` への像は `G` の正規部分群。
  - `exists_upperRamificationGroup_eq_lowerRamificationGroup_all`
    —— ★**すべての `m : ℝ`** で `G^m` は `ℕ` 添字の `G_n` に一致する
    (Y18 の版は `m ≥ 0` 限定だった。`m < 0` 側は `G^m = ⊤ = G_0` で埋まる)。
  - `normal_map_subtype_upperRamificationGroup`(上の 2 つの合成)。
  - `map_maximalIdeal_ne_bot_of_injective` / `faithfulSMul_subgroup` /
    `smulCommClass_subgroup_fixedRing`。
* §2 **具体層**:
  - ★★★`exists_sub_mem_maximalIdeal_inertiaFixedRing` ——
    **`L/L₀` は完全分岐**(`𝒪_L` の元はすべて `𝒪_{L₀} = 𝒪_L^{Gal(L/L₀)}` の元と
    `𝔪_L` を法として合同)。§1 の抽象核 + Y19b+c の `Gal(L/L₀) ≤ G_0`。
  - ★★★★`exists_uniformizer_adjoin_inertiaFixedRing_eq_top` ——
    **`𝒪_L = 𝒪_{L₀}[π]`**(Yoshida Lemma 5.11 を `L/L₀` に当てた形)。
    ★★**これが Y19 の言う「完全分岐の穴」の出口**である。
  - ★★★★`stage_upperRamification_coe_mul_coe_eq` ——
    **Yoshida Corollary 6.13(i) が段 `L/L₀` で使える**。
    ★★**Y18 `upperRamification_coe_mul_coe_eq` の仮定 10 個は全部供給できた**
    (Y19 は「できていない」と書いていた)。
  - `stageUpperRamification` —— 段の上付き分岐群 `Gal(L/L₀)^m ⊆ Gal(L/K)`。
    正規・反単調・`m ≤ 0` で `Gal(L/L₀)`・`m` 大で `⊥`。
* §4 **`Γ_K` への引き戻しと段データ**:
  - `StageGenerator` / `nonempty_stageGenerator`(開正規部分群は必ず `Gal(K̄/K(x))` の形)。
  - `absGalStage` —— 段データ `S N v`。`le_S` / `normal_S` / `antitone_S` が埋まった。
  - `absGalStageFiltration` / `ramificationFiltrationOfCompat` ——
    ★★**残る入力は `compat` ただ 1 つ**。

**★★埋まっていないもの(次のノード)**

★★**`compat`(`↑(S N' v) · N = ↑(S N v)`、`N' ≤ N`)は埋まっていない。**
★**近似(`trivialStageFiltration` / `inertiaStageFiltration`)で埋めることは
していない。** `ramificationFiltrationOfCompat` は `compat` を**仮定として受け取る**。

残っている数学を具体名で書く。`L ⊆ L′`、`G := Gal(L′/L′₀)`、
`H := G ⊓ Gal(L′/L) = Gal(L′/L·L′₀)` と置くと、本ファイルの
`stage_upperRamification_coe_mul_coe_eq` は
`G^{π′} m · H = G^{ϖ} m`(`ϖ` は `C := (𝒪_{L′})^H` の素元)を与える。
これを `compat` にするには、あと 2 つが要る:

1. ★**`C = (𝒪_{L′})^{Gal(L′/L·L′₀)} ≅ 𝒪_{L·L′₀}`** —— 固定環と中間体の整数環の同定。
   `Found/PGC/TotallyRamified.lean` の `adjoinIntegersRingHom` の像が
   ちょうど固定環であること。
2. ★★**不分岐底変換で上付き番号付けが変わらないこと** ——
   `Gal(L·L′₀/L′₀)^m` (`𝒪_{L·L′₀}` で計算)が `Gal(L/L₀)^m`(`𝒪_L` で計算)に
   対応すること。これは環同型 + 群同型に沿った `ramIndex` / `herbrandPhiGroup` の
   輸送であり、**本ファイルの外**である。
   ★Y19b+c が「合成体の Galois 群の同型 2 つは証明していない」と書いたものの、
   環レベルの相棒である。

★**Y19b+c が残したもう 1 つの穴(逆包含 `G_0 ≤ Gal(L/L₀)`)は本持ち場では要らなかった。**
使ったのは `Gal(L/L₀) ≤ G_0` の側だけである(§2)。

## ★退化の自己検査(近似と一致しないこと)

* `absGalStage_of_nonpos` : `v ≤ 0` で `S N v = I_K ⊔ N`。
* `exists_absGalStage_eq` : ★**`v` を大きくすると `S N v = N`**。
* `exists_absGalStage_ne_inertiaStageFiltration` : ★★したがって
  `I_K ≰ N`(= `L ⊄ K^ur`)なる `N` では、Y19b+c の
  `inertiaStageFiltration`(`v ≥ 0` で常に `I_K ⊔ N`)と**値が食い違う**。
  ⇒ 本ファイルの `absGalStage` は**上からの近似ではない**。
* Y19 の `trivialStageFiltration` は `v ≤ 0` で `⊤`、本ファイルは `I_K ⊔ N` なので、
  `I_K ⊔ N ≠ ⊤`(= `L/K` が完全分岐でない)なる `N` では**下からの退化とも異なる**。
* ★`v ≤ 0` の扱い: `Gal(L/L₀)^v = Gal(L/L₀)`(Y18 の `G^m = ⊤`、`m ≤ 0`)に合わせた。
  ★`v > 0` でも `S N v ⊇ N` なので、`Γ_K^v` は `v ≤ 0` で `I_K` を含む。
* ★`antitone` の向き: `Antitone` の定義どおり `v1 ≥ v2 → Γ^{v1} ⊆ Γ^{v2}`。
  ★上付きは番号が大きいほど小さい。
* ★`ℕ∞` の切り詰め引き算・除算は 1 箇所も書いていない(`lean-idioms.md` #102)。
* ★中間体の塔(`↥L` の上の `↥L′`)は 1 つも作っていない(#59 回避。Y19b+c と同じ流儀)。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. ★**一意化元 `π` を選択公理で 1 つ選んでいる**(`stageUniformizer`)。
   `G^m` は `π` の取り方に依らない(`G_n` が `π` に依らないから)が、
   **その独立性は本ファイルでは証明していない**。段データが `π` の選択に依存して
   見えるのはそのためである。★消費側(`limit` / `compat`)には影響しない。
2. ★同様に、開正規部分群 `N` に対する原始元 `x`(`StageGenerator`)も選択している。
   異なる `x` が同じ `S N v` を与えることは**証明していない**。
3. ★原文の `v` の範囲は `v > 0` だが、`Interface` の `RamificationFiltration` に
   合わせて全実数で定義する(Y19 の逸脱 1 を引き継いだだけ)。
4. ★本ファイルは**原典が名前を付けていない組み立てノード**なので `.src` を持たない
   (`InertiaReduction.lean` と同じ理由)。使っている原典項目は
   Yoshida Definition 6.12 / Corollary 6.13(i) / Lemma 5.11 であり、
   `.src` はそれぞれ `UpperRamificationGroup.lean` / `FixedRingMonogenic.lean` にある。
5. ★**既存の `Found/PGC/*.lean` は 1 行も書き換えていない**(import のみ)。
-/

namespace ABC3.Found.PGC

open scoped Pointwise

/-! ## §1 抽象核 —— 純可換環 + 群作用(分岐・付値・Galois の語彙は 1 つも出てこない) -/

open IsLocalRing in
/-- ★★★**抽象核** —— `H` が剰余体に自明に作用するなら、`B` のどの元も
**`H` 不変な元と `𝔪` を法として合同**である。

証明は **Teichmüller 代表元の一意性**: `b` の剰余が `0` でなければ
`ζ^{q-1} = 1` かつ `ζ ≡ b` なる `ζ` が一意に取れる(`ABC3.Found.exists_teichmullerRep` /
`teichmullerRep_unique`)。`σ ∈ H` に対して `σ ζ` も同じ条件を満たすので `σ ζ = ζ`。

★★これが原典 §6.1 の "totally ramified"(Y18 の `hres` / `hresA`)の中身である。
★**分岐・付値・Galois の語彙が 1 つも出てこない**(可換環 `B` と作用する群 `G` だけ)。
★`HenselianLocalRing` と剰余体の有限性は落とせない(Teichmüller が要る)。 -/
theorem exists_fixed_sub_mem_maximalIdeal {B : Type*} [CommRing B] [HenselianLocalRing B]
    [Fintype (IsLocalRing.ResidueField B)] {G : Type*} [Group G] [MulSemiringAction G B]
    {H : Subgroup G} (h : ∀ σ ∈ H, ∀ y : B, σ • y - y ∈ maximalIdeal B) (b : B) :
    ∃ c : B, (∀ σ ∈ H, σ • c = c) ∧ b - c ∈ maximalIdeal B := by
  by_cases hb : IsLocalRing.residue B b = 0
  · refine ⟨0, fun σ _ => smul_zero σ, ?_⟩
    rw [sub_zero, ← IsLocalRing.residue_eq_zero_iff]
    exact hb
  · obtain ⟨ζ, hζpow, hζres⟩ := ABC3.Found.exists_teichmullerRep (IsLocalRing.residue B b) hb
    refine ⟨ζ, fun σ hσ => ?_, ?_⟩
    · refine ABC3.Found.teichmullerRep_unique ?_ hζpow ?_
      · rw [← smul_pow', hζpow, smul_one]
      · rw [← sub_eq_zero, ← map_sub, IsLocalRing.residue_eq_zero_iff]
        exact h σ hσ ζ
    · rw [← IsLocalRing.residue_eq_zero_iff, map_sub, hζres, sub_self]

open IsLocalRing in
/-- ★★**上の固定環版** —— `H ≤ 𝔪.inertia G`(= `H ⊆ G_0`)なら
`B^H` の剰余体は `B` の剰余体に等しい。★これがそのまま Y18 の `hres` の形である。 -/
theorem exists_sub_mem_maximalIdeal_fixedRing_of_le_inertia {B : Type*} [CommRing B]
    [HenselianLocalRing B] [Fintype (IsLocalRing.ResidueField B)]
    {G : Type*} [Group G] [MulSemiringAction G B] {H : Subgroup G}
    (h : H ≤ (maximalIdeal B).inertia G) (b : B) :
    ∃ c : ↥(fixedRing B H), b - algebraMap (↥(fixedRing B H)) B c ∈ maximalIdeal B := by
  obtain ⟨c, hc, hb⟩ := exists_fixed_sub_mem_maximalIdeal
    (H := H) (fun σ hσ y => AddSubgroup.mem_inertia.mp (h hσ) y) b
  exact ⟨⟨c, mem_fixedRing.2 hc⟩, hb⟩

open IsLocalRing in
/-- ★**抽象核** —— 群作用は `𝔪^k` を保つ(環自己同型は極大イデアルを保つから)。 -/
theorem smul_mem_pow_maximalIdeal {B : Type*} [CommRing B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] (σ : G) {k : ℕ} {x : B}
    (hx : x ∈ (maximalIdeal B) ^ k) : σ • x ∈ (maximalIdeal B) ^ k := by
  have hle : Ideal.map (MulSemiringAction.toRingHom G B σ) (maximalIdeal B) ≤ maximalIdeal B := by
    rw [Ideal.map_le_iff_le_comap]
    intro y hy
    exact smul_mem_maximalIdeal σ hy
  have h1 : σ • x ∈ Ideal.map (MulSemiringAction.toRingHom G B σ) ((maximalIdeal B) ^ k) :=
    Ideal.mem_map_of_mem _ hx
  rw [Ideal.map_pow] at h1
  exact Ideal.pow_right_mono hle k h1

/-- ★★**抽象核** —— `H ⊴ G` なら `G_n(H)` の `G` への像は `G` の**正規**部分群。

`τ σ τ⁻¹` が `x` に及ぼす差は `τ • (σ • (τ⁻¹ • x) − τ⁻¹ • x)` であり、
`𝔪^{n+1}` は `τ` で保たれる。★分岐も付値も出てこない。 -/
theorem normal_map_subtype_lowerRamificationGroup {B : Type*} [CommRing B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {H : Subgroup G} [hH : H.Normal] (n : ℕ) :
    (Subgroup.map H.subtype (lowerRamificationGroup B (H : Type _) n)).Normal := by
  constructor
  rintro _ ⟨⟨s, hs⟩, hmem, rfl⟩ τ
  simp only [SetLike.mem_coe, mem_lowerRamificationGroup_iff_forall] at hmem
  refine ⟨⟨τ * s * τ⁻¹, hH.conj_mem s hs τ⟩, ?_, rfl⟩
  rw [SetLike.mem_coe, mem_lowerRamificationGroup_iff_forall]
  intro y
  show (τ * s * τ⁻¹) • y - y ∈ _
  have hrw : (τ * s * τ⁻¹) • y - y = τ • (s • (τ⁻¹ • y) - τ⁻¹ • y) := by
    rw [smul_sub, ← mul_smul, ← mul_smul, ← mul_smul]
    simp
  rw [hrw]
  exact smul_mem_pow_maximalIdeal τ (hmem _)

open IsLocalRing in
/-- `r ≤ 0` では実数添字の下付き分岐群は `⊤`(`G_0 = G` と反単調性)。 -/
theorem ramificationGroupReal_eq_top_of_nonpos {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] {α : B}
    (huni : maximalIdeal B = Ideal.span {α}) {r : ℝ} (hr : r ≤ 0) :
    ramificationGroupReal (G := G) α r = ⊤ := by
  refine eq_top_iff.2 ?_
  rw [← ramificationGroupReal_zero_eq_top (G := G) huni]
  exact ramificationGroupReal_antitone α hr

open IsLocalRing in
/-- ★★**すべての `m : ℝ`** で `G^m` は `ℕ` 添字の下付き分岐群に一致する。

Y18 の `exists_upperRamificationGroup_eq_lowerRamificationGroup` は `m ≥ 0` 限定だった。
`ψ_G(m) < 0` の側は `G^m = ⊤ = G_0` で埋まる(`ramificationGroupReal_eq_top_of_nonpos` +
`lowerRamificationGroup_zero_eq_top`)。★これが `G^m` の正規性を出すための入口である。 -/
theorem exists_upperRamificationGroup_eq_lowerRamificationGroup_all {A B : Type*} [CommRing A]
    [CommRing B] [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B] {α : B}
    (huni : maximalIdeal B = Ideal.span {α}) (hadj : Algebra.adjoin A ({α} : Set B) = ⊤)
    (m : ℝ) : ∃ n : ℕ, upperRamificationGroup G α m = lowerRamificationGroup B G n := by
  by_cases hm : 0 ≤ herbrandPsiGroup G α m
  · refine ⟨⌈herbrandPsiGroup G α m⌉₊,
      upperRamificationGroup_eq_lowerRamificationGroup (A := A) huni hadj _ ?_ (Nat.le_ceil _)⟩
    have h := Nat.ceil_lt_add_one hm
    linarith
  · replace hm : herbrandPsiGroup G α m < 0 := not_le.mp hm
    refine ⟨0, ?_⟩
    rw [upperRamificationGroup_def, ramificationGroupReal_eq_top_of_nonpos huni hm.le,
      lowerRamificationGroup_zero_eq_top (A := A)
        (by rw [huni]; exact Ideal.mem_span_singleton_self α) hadj]

open IsLocalRing in
/-- ★★**`H ⊴ G` なら `H^m` の `G` への像は `G` の正規部分群**(上の 2 本の合成)。 -/
theorem normal_map_subtype_upperRamificationGroup {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {H : Subgroup G} [H.Normal] [Fintype H]
    [SMulCommClass (H : Type _) A B] {α : B}
    (huni : maximalIdeal B = Ideal.span {α}) (hadj : Algebra.adjoin A ({α} : Set B) = ⊤)
    (m : ℝ) :
    (Subgroup.map H.subtype (upperRamificationGroup (H : Type _) α m)).Normal := by
  obtain ⟨n, hn⟩ :=
    exists_upperRamificationGroup_eq_lowerRamificationGroup_all (A := A) (G := (H : Type _))
      huni hadj m
  rw [hn]
  exact normal_map_subtype_lowerRamificationGroup n

open IsLocalRing in
/-- `A` が DVR で `A → B` が単射なら `𝔪_A` の像は `0` でない
(`exists_uniformizer_adjoin_eq` の `hne` を供給する)。 -/
theorem map_maximalIdeal_ne_bot_of_injective {A B : Type*} [CommRing A] [IsDomain A]
    [IsDiscreteValuationRing A] [CommRing B] [Algebra A B]
    (hinj : Function.Injective (algebraMap A B)) :
    Ideal.map (algebraMap A B) (maximalIdeal A) ≠ ⊥ := by
  obtain ⟨a, ha, ha0⟩ := exists_mem_maximalIdeal_map_ne_zero (B := B) hinj
  intro h
  refine ha0 ?_
  have hmem : algebraMap A B a ∈ Ideal.map (algebraMap A B) (maximalIdeal A) :=
    Ideal.mem_map_of_mem _ ha
  rw [h, Ideal.mem_bot] at hmem
  exact hmem

/-- 部分群への制限も忠実。 -/
instance faithfulSMul_subgroup {α G : Type*} [Group G] [MulAction G α] [FaithfulSMul G α]
    (H : Subgroup G) : FaithfulSMul (H : Type _) α :=
  ⟨fun h => Subtype.ext (eq_of_smul_eq_smul h)⟩

/-- `H` の元は固定環の元と可換に作用する。 -/
instance smulCommClass_subgroup_fixedRing {B : Type*} [CommRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] (H : Subgroup G) :
    SMulCommClass (H : Type _) ↥(fixedRing B H) B where
  smul_comm σ a b := by
    show (σ : G) • ((a : B) * b) = (a : B) * ((σ : G) • b)
    rw [smul_mul', mem_fixedRing.1 a.2 (σ : G) σ.2]

/-! ## §2 具体層 -/

open ABC3.Skeleton.PGC ABC3.Interface.PGC
open scoped NNReal Valued

attribute [local instance] isDiscreteValuationRing_adjoinIntegers

variable {p : ℕ} [Fact p.Prime]

/-- `Gal(L/L₀) = I_K の像` は `Gal(L/K)` の正規部分群(`I_K ⊴ Γ_K` の像だから)。
★`lean-idioms.md` #126 の教訓により、証明の中の `haveI` では間に合わないので
**前もって instance にしておく**。 -/
instance normal_inertiaGal (K : PAdicLocalField p)
    (L : IntermediateField K.carrier K.closure) [Normal K.carrier L] :
    (inertiaGal K L).Normal := by
  haveI := isGalois_closure K
  exact (normal_absInertia K).map _ (surjective_restrictNormalHom L)

/-- `L = K(x)` の惰性部分群 `Gal(L/L₀)`。★`abbrev` なのでインスタンス探索が透過する。 -/
noncomputable abbrev inertiaGalAdjoin (K : PAdicLocalField p) (x : K.closure)
    [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    Subgroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :=
  inertiaGal K (IntermediateField.adjoin K.carrier ({x} : Set K.closure))

section Stage

variable (K : PAdicLocalField p) (x : K.closure)
  [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
  [Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]

/-- `Gal(L/K)` は有限群(`L/K` が有限次だから)。 -/
noncomputable instance fintypeGalAdjoin :
    Fintype ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :=
  Fintype.ofFinite _

/-- 惰性群も有限群。 -/
noncomputable instance fintypeInertiaGalAdjoin : Fintype ↥(inertiaGalAdjoin K x) :=
  Fintype.ofFinite _

/-- `𝒪_L` は `𝒪_{L₀} = 𝒪_L^{Gal(L/L₀)}` 上有限。 -/
theorem module_finite_inertiaFixedRing :
    Module.Finite ↥(fixedRing (adjoinIntegers K x) (inertiaGalAdjoin K x))
      (adjoinIntegers K x) := by
  haveI := module_finite_adjoinIntegers K x
  letI := fixedRingAlgebra 𝒪[K.carrier] (adjoinIntegers K x) (inertiaGalAdjoin K x)
  haveI : IsScalarTower 𝒪[K.carrier] ↥(fixedRing (adjoinIntegers K x) (inertiaGalAdjoin K x))
      (adjoinIntegers K x) := isScalarTower_fixedRing
  exact Module.Finite.of_restrictScalars_finite 𝒪[K.carrier] _ _

/-- ★★★**`L/L₀` は完全分岐** —— 剰余体が伸びない形。

`𝒪_L` のどの元も `𝒪_{L₀} = 𝒪_L^{Gal(L/L₀)}` の元と `𝔪_L` を法として合同である。
★これが Y18 の 10 個の仮定のうち `hres` / `hresA` の中身であり、
Y19 が「供給できていない」と書いたものである。

証明は §1 の抽象核 `exists_sub_mem_maximalIdeal_fixedRing_of_le_inertia`
(Teichmüller 代表元の一意性)に、Y19b+c の
`inertiaGal_le_lowerRamificationGroupAdjoin_zero`(`Gal(L/L₀) ≤ G_0`)を代入するだけ。
★★**Y19b+c が残した「逆包含 `G_0 ≤ Gal(L/L₀)`」は要らない。** -/
theorem exists_sub_mem_maximalIdeal_inertiaFixedRing (b : adjoinIntegers K x) :
    ∃ c : ↥(fixedRing (adjoinIntegers K x) (inertiaGalAdjoin K x)),
      b - algebraMap ↥(fixedRing (adjoinIntegers K x) (inertiaGalAdjoin K x))
        (adjoinIntegers K x) c ∈ IsLocalRing.maximalIdeal (adjoinIntegers K x) := by
  haveI : Fintype (IsLocalRing.ResidueField (adjoinIntegers K x)) := Fintype.ofFinite _
  refine exists_sub_mem_maximalIdeal_fixedRing_of_le_inertia ?_ b
  rw [← lowerRamificationGroup_zero_eq_inertia]
  exact inertiaGal_le_lowerRamificationGroupAdjoin_zero K x

/-- ★★★★**Yoshida Lemma 5.11 を `L/L₀` に当てた形** —— `𝒪_L = 𝒪_{L₀}[π]`。

★★これが Y19 の言う「完全分岐の穴」の出口である: 一般の有限 Galois `L/K` に対して
底を `L₀ = L ⊓ K^ur` へ移せば、Yoshida §6.1 の設定(`Algebra.adjoin A {π} = ⊤`)が
**実際に満たされる**。 -/
theorem exists_uniformizer_adjoin_inertiaFixedRing_eq_top :
    ∃ π : adjoinIntegers K x,
      IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {π} ∧
        Algebra.adjoin ↥(fixedRing (adjoinIntegers K x) (inertiaGalAdjoin K x))
          ({π} : Set (adjoinIntegers K x)) = ⊤ := by
  haveI := isDiscreteValuationRing_adjoinIntegers K x
  haveI : Fintype ↥(inertiaGalAdjoin K x) := Fintype.ofFinite _
  haveI := module_finite_inertiaFixedRing K x
  exact exists_uniformizer_adjoin_eq (exists_sub_mem_maximalIdeal_inertiaFixedRing K x)
    (map_maximalIdeal_ne_bot_of_injective fixedRing_injective)

/-- 段の一意化元(選択)。★`𝔪_L = (π)` かつ `𝒪_L = 𝒪_{L₀}[π]`。 -/
noncomputable def stageUniformizer : adjoinIntegers K x :=
  (exists_uniformizer_adjoin_inertiaFixedRing_eq_top K x).choose

theorem maximalIdeal_eq_span_stageUniformizer :
    IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {stageUniformizer K x} :=
  (exists_uniformizer_adjoin_inertiaFixedRing_eq_top K x).choose_spec.1

theorem adjoin_stageUniformizer_eq_top :
    Algebra.adjoin ↥(fixedRing (adjoinIntegers K x) (inertiaGalAdjoin K x))
      ({stageUniformizer K x} : Set (adjoinIntegers K x)) = ⊤ :=
  (exists_uniformizer_adjoin_inertiaFixedRing_eq_top K x).choose_spec.2

/-- `𝒪_L` は `𝒪_{L₀}` 上ネーター加群。 -/
theorem isNoetherian_inertiaFixedRing :
    IsNoetherian ↥(fixedRing (adjoinIntegers K x) (inertiaGalAdjoin K x)) (adjoinIntegers K x) := by
  haveI := isDiscreteValuationRing_adjoinIntegers K x
  haveI := module_finite_inertiaFixedRing K x
  exact isNoetherian_of_isNoetherianRing_of_finite _ _

/-- ★★★★**Yoshida Corollary 6.13(i) を段 `L/L₀` に当てた形**。

★★これが本ファイルの**答え**である: Y18 `upperRamification_coe_mul_coe_eq` の
**10 個の仮定はすべて供給できた**。決め手は §2 の
`adjoin_stageUniformizer_eq_top`(`𝒪_L = 𝒪_{L₀}[π]`)であり、それは
Y19b+c の `Gal(L/L₀) ≤ G_0` と §1 の Teichmüller 抽象核から出た。

* `hcomp` = `algebraMap_smul_fixedRing`、`hHtriv` = `smul_fixedRing_eq_self`、
* `hinj` = `fixedRing_injective`、`hfixC` = `exists_algebraMap_fixedRing`、
* `hres` = `exists_sub_mem_fixedRing` + §2 の完全分岐、
* `hAC` = `exists_algebraMap_fixedRing_eq`、`hϖ` = `rfl`、
* `hadj` = §2 の `adjoin_stageUniformizer_eq_top`(★**Y19 が塞がれていると書いたもの**)、
* `hfix` = Y20 の `fixedRing_mem_adjoin_uniformizer` + §2 の完全分岐。 -/
theorem stage_upperRamification_coe_mul_coe_eq
    (H : Subgroup ↥(inertiaGalAdjoin K x)) [H.Normal] [Fintype ↥H]
    {ϖ : ↥(fixedRing (adjoinIntegers K x) H)} (hϖ : Irreducible ϖ) (m : ℝ) :
    ((upperRamificationGroup ↥(inertiaGalAdjoin K x) (stageUniformizer K x) m :
        Subgroup ↥(inertiaGalAdjoin K x)) : Set ↥(inertiaGalAdjoin K x)) * (H : Set _)
      = ((upperRamificationGroup ↥(inertiaGalAdjoin K x) ϖ m :
          Subgroup ↥(inertiaGalAdjoin K x)) : Set ↥(inertiaGalAdjoin K x)) := by
  haveI := isDiscreteValuationRing_adjoinIntegers K x
  haveI := module_finite_inertiaFixedRing K x
  haveI := isNoetherian_inertiaFixedRing K x
  have hresA := exists_sub_mem_maximalIdeal_inertiaFixedRing K x
  have hAne : ∃ a ∈ IsLocalRing.maximalIdeal
      ↥(fixedRing (adjoinIntegers K x) (inertiaGalAdjoin K x)),
      algebraMap ↥(fixedRing (adjoinIntegers K x) (inertiaGalAdjoin K x))
        (adjoinIntegers K x) a ≠ 0 :=
    exists_mem_maximalIdeal_map_ne_zero fixedRing_injective
  exact upperRamification_coe_mul_coe_eq
    (A := ↥(fixedRing (adjoinIntegers K x) (inertiaGalAdjoin K x)))
    (fun ρ c => algebraMap_smul_fixedRing ρ c) smul_fixedRing_eq_self
    ((IsDiscreteValuationRing.irreducible_iff_uniformizer _).2
      (maximalIdeal_eq_span_stageUniformizer K x))
    fixedRing_injective exists_algebraMap_fixedRing
    (fun b => exists_sub_mem_fixedRing hresA b) exists_algebraMap_fixedRing_eq rfl
    (adjoin_stageUniformizer_eq_top K x)
    (fixedRing_mem_adjoin_uniformizer hresA hAne hϖ) m

/-! ### §3 段の上付き分岐群 -/

/-- ★★★★**段の上付き分岐群** `Gal(L/L₀)^m ⊆ Gal(L/K)`。

Yoshida Definition 6.12 を、完全分岐拡大 `L/L₀`(§2 で設定が満たされることを示した)に
当てて `Gal(L/K)` の中へ押し出したもの。

★★**近似ではない**: `m ≤ 0` では `Gal(L/L₀)`(`stageUpperRamification_of_nonpos`)、
`m` が大きいところでは `⊥`(`exists_stageUpperRamification_eq_bot`)。
Y19b+c の `inertiaStageFiltration` は `v ≥ 0` で常に `I_K` なので**別物**であり、
Y19 の `trivialStageFiltration` は `v > 0` で常に `⊥` なのでこれも**別物**である。 -/
noncomputable def stageUpperRamification (m : ℝ) :
    Subgroup ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :=
  Subgroup.map (inertiaGalAdjoin K x).subtype
    (upperRamificationGroup (↥(inertiaGalAdjoin K x)) (stageUniformizer K x) m)

/-- 段の分岐群は惰性群に含まれる。 -/
theorem stageUpperRamification_le_inertiaGal (m : ℝ) :
    stageUpperRamification K x m ≤ inertiaGalAdjoin K x :=
  Subgroup.map_subtype_le _

/-- ★段の分岐群は `Gal(L/K)` の**正規**部分群。 -/
instance normal_stageUpperRamification (m : ℝ) : (stageUpperRamification K x m).Normal :=
  normal_map_subtype_upperRamificationGroup
    (A := ↥(fixedRing (adjoinIntegers K x) (inertiaGalAdjoin K x)))
    (maximalIdeal_eq_span_stageUniformizer K x) (adjoin_stageUniformizer_eq_top K x) m

/-- ★段の分岐群は `m` について**反単調**(上付きは番号が大きいほど小さい)。 -/
theorem antitone_stageUpperRamification : Antitone (stageUpperRamification K x) := fun _ _ hab =>
  Subgroup.map_mono
    (upperRamificationGroup_antitone (maximalIdeal_eq_span_stageUniformizer K x) hab)

/-- ★**`m ≤ 0` では惰性群そのもの**(Y18 の `G^m = ⊤` を押し出した形)。 -/
theorem stageUpperRamification_of_nonpos {m : ℝ} (hm : m ≤ 0) :
    stageUpperRamification K x m = inertiaGalAdjoin K x := by
  have htop : upperRamificationGroup (↥(inertiaGalAdjoin K x)) (stageUniformizer K x) m = ⊤ := by
    refine eq_top_iff.2 ?_
    rw [← upperRamificationGroup_zero (maximalIdeal_eq_span_stageUniformizer K x)]
    exact upperRamificationGroup_antitone (maximalIdeal_eq_span_stageUniformizer K x) hm
  rw [stageUpperRamification, htop]
  refine le_antisymm (Subgroup.map_subtype_le _) ?_
  intro g hg
  exact ⟨⟨g, hg⟩, Subgroup.mem_top _, rfl⟩

/-- ★★**退化していないことの検査(上から)** —— `m` を十分大きく取れば `⊥` になる。

★これで `Y19b+c` の `inertiaStageFiltration`(`v ≥ 0` で常に `I_K`)と**一致しない**
ことが分かる(惰性群が自明でない `L` を取れば `m = 0` と大きい `m` で値が違う)。 -/
theorem exists_stageUpperRamification_eq_bot : ∃ m : ℝ, stageUpperRamification K x m = ⊥ := by
  haveI : IsNoetherianRing (adjoinIntegers K x) := inferInstance
  obtain ⟨N, hN⟩ := exists_lowerRamificationGroup_eq_bot
    (A := ↥(fixedRing (adjoinIntegers K x) (inertiaGalAdjoin K x)))
    (G := ↥(inertiaGalAdjoin K x)) (adjoin_stageUniformizer_eq_top K x)
  refine ⟨herbrandPhiGroup (↥(inertiaGalAdjoin K x)) (stageUniformizer K x) (N : ℝ), ?_⟩
  rw [stageUpperRamification, upperRamificationGroup_herbrandPhiGroup,
    ramificationGroupReal_eq_of_mem_Ioc
      (A := ↥(fixedRing (adjoinIntegers K x) (inertiaGalAdjoin K x)))
      (maximalIdeal_eq_span_stageUniformizer K x) (adjoin_stageUniformizer_eq_top K x) N
      (by linarith) le_rfl,
    hN N le_rfl, Subgroup.map_bot]

end Stage

/-! ## §4 `Γ_K` への引き戻しと段データ -/

/-- 開正規部分群 `N ⊴ Γ_K` の「段の生成元」——
`N = Gal(K̄/K(x))` かつ `K(x)/K` は有限次正規。

★`FiniteDimensional` / `Normal` は `stageUpperRamification` の**型**に現れるので、
`∃` 文ではなく構造体で持ち回る(`lean-idioms.md` #126)。 -/
structure StageGenerator (K : PAdicLocalField p) (N : Subgroup K.absGal) where
  /-- 原始元 -/
  gen : K.closure
  /-- `K(x)/K` は有限次 -/
  finiteDimensional :
    FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({gen} : Set K.closure))
  /-- `K(x)/K` は正規 -/
  normal : Normal K.carrier (IntermediateField.adjoin K.carrier ({gen} : Set K.closure))
  /-- `N` は `K(x)` の固定化部分群 -/
  fixingSubgroup_eq :
    (IntermediateField.adjoin K.carrier ({gen} : Set K.closure)).fixingSubgroup = N

/-- ★**開正規部分群には必ず段の生成元がある**(原始元定理 + 無限次 Galois 対応)。 -/
theorem nonempty_stageGenerator (K : PAdicLocalField p) {N : Subgroup K.absGal}
    (hN : N ∈ openNormalBase K.absGal) : Nonempty (StageGenerator K N) := by
  haveI := isGalois_closure K
  haveI := hN.2
  obtain ⟨x, hx⟩ := exists_adjoin_eq_fixedField K N hN.1
  have hclosed : IsClosed (N : Set K.absGal) := Subgroup.isClosed_of_isOpen N hN.1
  have hfix : (IntermediateField.fixedField N).fixingSubgroup = N :=
    InfiniteGalois.fixingSubgroup_fixedField (⟨N, hclosed⟩ : ClosedSubgroup K.absGal)
  have hfin : FiniteDimensional K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
    rw [hx]; exact finiteDimensional_fixedField_of_isOpen K N hN.1
  have hgal : IsGalois K.carrier (IntermediateField.fixedField N) :=
    (InfiniteGalois.normal_iff_isGalois _).mp (by rw [hfix]; exact hN.2)
  have hnor : Normal K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
    haveI := hgal
    rw [hx]
    infer_instance
  exact ⟨⟨x, hfin, hnor, by rw [hx, hfix]⟩⟩

/-- 段の生成元から作る `Γ_K` の部分群 —— `Gal(L/L₀)^v` の `Γ_K` への引き戻し。 -/
noncomputable def StageGenerator.stage {K : PAdicLocalField p} {N : Subgroup K.absGal}
    (g : StageGenerator K N) (v : ℝ) : Subgroup K.absGal :=
  letI := g.finiteDimensional
  letI := g.normal
  Subgroup.comap (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      ((IntermediateField.adjoin K.carrier ({g.gen} : Set K.closure)) : Type _))
    (stageUpperRamification K g.gen v)

/-- ★★**段データの本体** `S N v` —— `N` が開正規なら `Gal(L/L₀)^v` の引き戻し、
そうでなければ `⊤`(`base` の外なので `limit` には効かない)。 -/
noncomputable def absGalStage (K : PAdicLocalField p) (N : Subgroup K.absGal) (v : ℝ) :
    Subgroup K.absGal :=
  letI := Classical.dec (Nonempty (StageGenerator K N))
  if h : Nonempty (StageGenerator K N) then h.some.stage v else ⊤

theorem absGalStage_eq_stage {K : PAdicLocalField p} {N : Subgroup K.absGal}
    (h : Nonempty (StageGenerator K N)) (v : ℝ) : absGalStage K N v = h.some.stage v := by
  rw [absGalStage]
  exact dif_pos h

/-- `N ≤ S N v`(引き戻しなので `ker` を含む)。 -/
theorem le_stage {K : PAdicLocalField p} {N : Subgroup K.absGal} (g : StageGenerator K N)
    (v : ℝ) : N ≤ g.stage v := by
  letI := g.finiteDimensional
  letI := g.normal
  intro σ hσ
  have h1 : σ ∈ (IntermediateField.adjoin K.carrier ({g.gen} : Set K.closure)).fixingSubgroup := by
    rw [g.fixingSubgroup_eq]; exact hσ
  rw [← ker_restrictNormalHom_eq_fixingSubgroup
    (E := K.closure) (IntermediateField.adjoin K.carrier ({g.gen} : Set K.closure)),
    MonoidHom.mem_ker] at h1
  show σ ∈ Subgroup.comap _ _
  rw [Subgroup.mem_comap, h1]
  exact one_mem _

theorem normal_stage {K : PAdicLocalField p} {N : Subgroup K.absGal} (g : StageGenerator K N)
    (v : ℝ) : (g.stage v).Normal := by
  letI := g.finiteDimensional
  letI := g.normal
  exact Subgroup.Normal.comap (normal_stageUpperRamification K g.gen v) _

theorem antitone_stage {K : PAdicLocalField p} {N : Subgroup K.absGal} (g : StageGenerator K N) :
    Antitone g.stage := by
  letI := g.finiteDimensional
  letI := g.normal
  exact fun _ _ hab => Subgroup.comap_mono (antitone_stageUpperRamification K g.gen hab)

theorem le_absGalStage (K : PAdicLocalField p) {N : Subgroup K.absGal}
    (hN : N ∈ openNormalBase K.absGal) (v : ℝ) : N ≤ absGalStage K N v := by
  rw [absGalStage_eq_stage (nonempty_stageGenerator K hN)]
  exact le_stage _ v

theorem normal_absGalStage (K : PAdicLocalField p) {N : Subgroup K.absGal}
    (hN : N ∈ openNormalBase K.absGal) (v : ℝ) : (absGalStage K N v).Normal := by
  rw [absGalStage_eq_stage (nonempty_stageGenerator K hN)]
  exact normal_stage _ v

theorem antitone_absGalStage (K : PAdicLocalField p) {N : Subgroup K.absGal}
    (hN : N ∈ openNormalBase K.absGal) : Antitone (absGalStage K N) := by
  intro a b hab
  rw [absGalStage_eq_stage (nonempty_stageGenerator K hN),
    absGalStage_eq_stage (nonempty_stageGenerator K hN)]
  exact antitone_stage _ hab

/-! ### 退化の自己検査 —— 近似の 2 つと一致しないこと -/

/-- ★★**`v ≤ 0` では `I_K ⊔ N`** —— Y19b+c の `comap_inertiaGal` がそのまま効く。 -/
theorem stage_of_nonpos {K : PAdicLocalField p} {N : Subgroup K.absGal} (g : StageGenerator K N)
    {v : ℝ} (hv : v ≤ 0) : g.stage v = absInertia K ⊔ N := by
  letI := g.finiteDimensional
  letI := g.normal
  rw [StageGenerator.stage, stageUpperRamification_of_nonpos K g.gen hv]
  rw [show (inertiaGalAdjoin K g.gen) =
      inertiaGal K (IntermediateField.adjoin K.carrier ({g.gen} : Set K.closure)) from rfl,
    comap_inertiaGal, g.fixingSubgroup_eq]

/-- ★★**`v` を大きくすると `N` そのもの** —— `Gal(L/L₀)^v = ⊥` になるから。

★★これが「近似ではない」ことの決め手である: Y19b+c の `inertiaStageFiltration` は
`v ≥ 0` で常に `I_K ⊔ N` なので、`I_K ≰ N`(= `L ⊄ K^ur`)なる `N` では
**値が食い違う**。 -/
theorem exists_stage_eq {K : PAdicLocalField p} {N : Subgroup K.absGal}
    (g : StageGenerator K N) : ∃ v : ℝ, g.stage v = N := by
  letI := g.finiteDimensional
  letI := g.normal
  obtain ⟨v, hv⟩ := exists_stageUpperRamification_eq_bot K g.gen
  refine ⟨v, ?_⟩
  rw [StageGenerator.stage, hv, MonoidHom.comap_bot, ker_restrictNormalHom_eq_fixingSubgroup,
    g.fixingSubgroup_eq]

/-- ★`v ≤ 0` での値。 -/
theorem absGalStage_of_nonpos (K : PAdicLocalField p) {N : Subgroup K.absGal}
    (hN : N ∈ openNormalBase K.absGal) {v : ℝ} (hv : v ≤ 0) :
    absGalStage K N v = absInertia K ⊔ N := by
  rw [absGalStage_eq_stage (nonempty_stageGenerator K hN)]
  exact stage_of_nonpos _ hv

/-- ★大きい `v` での値。 -/
theorem exists_absGalStage_eq (K : PAdicLocalField p) {N : Subgroup K.absGal}
    (hN : N ∈ openNormalBase K.absGal) : ∃ v : ℝ, absGalStage K N v = N := by
  obtain ⟨v, hv⟩ := exists_stage_eq (nonempty_stageGenerator K hN).some
  exact ⟨v, by rw [absGalStage_eq_stage (nonempty_stageGenerator K hN)]; exact hv⟩

/-- ★★★**`inertiaStageFiltration`(上からの近似)と一致しない**。

`I_K ≰ N` なら、`v` を大きく取ったところで
`absGalStage K N v = N ≠ I_K ⊔ N = (inertiaStageFiltration K).S N v` である。 -/
theorem exists_absGalStage_ne_inertiaStageFiltration (K : PAdicLocalField p)
    {N : Subgroup K.absGal} (hN : N ∈ openNormalBase K.absGal) (hne : ¬ absInertia K ≤ N) :
    ∃ v : ℝ, absGalStage K N v ≠ (inertiaStageFiltration K).S N v := by
  obtain ⟨v, hv⟩ := exists_absGalStage_eq K hN
  refine ⟨max v 0, ?_⟩
  have hv0 : ¬ (max v 0 : ℝ) < 0 := by simp
  have hmax : absGalStage K N (max v 0) = N := by
    refine le_antisymm ?_ (le_absGalStage K hN _)
    calc absGalStage K N (max v 0)
        ≤ absGalStage K N v := antitone_absGalStage K hN (le_max_left v 0)
      _ = N := hv
  rw [hmax]
  show N ≠ (if (max v 0 : ℝ) < 0 then (⊤ : Subgroup K.absGal) else absInertia K ⊔ N)
  simp only [hv0, if_false]
  intro hcon
  exact hne (le_trans le_sup_left (le_of_eq hcon.symm))

/-! ### ★★段データの組み立て —— 残るのは `compat` だけ -/

/-- ★★★**`compat` を仮定として受け取る段データ**。

`le_S` / `normal_S` / `antitone_S` は本ファイルで埋まっている。
★**残っているのは `compat`(Yoshida Corollary 6.13(i) を有限段の対に当てる部分)だけ**
であり、それが次のノードである。 -/
noncomputable def absGalStageFiltration (K : PAdicLocalField p)
    (hcompat : ∀ ⦃M⦄, M ∈ openNormalBase K.absGal → ∀ ⦃N⦄, N ∈ openNormalBase K.absGal →
      N ≤ M → ∀ v, ((absGalStage K N v : Subgroup K.absGal) : Set K.absGal)
        * (M : Set K.absGal) = ((absGalStage K M v : Subgroup K.absGal) : Set K.absGal)) :
    StageFiltration K.absGal :=
  StageFiltration.ofOpenNormal (absGalStage K) (fun _ hN v => le_absGalStage K hN v)
    (fun _ hN v => normal_absGalStage K hN v) (fun _ hN => antitone_absGalStage K hN) hcompat

@[simp] theorem absGalStageFiltration_S (K : PAdicLocalField p) (hcompat) :
    (absGalStageFiltration K hcompat).S = absGalStage K := rfl

/-- ★★★**`RamificationFiltration p` の本物** —— `compat` を仮定として受け取る形。

★★`Skeleton/PGC/Section2.lean` の `prop_2_1` / `prop_2_2` はこれを入力に取れる。
★★**`compat` が埋まるまでは仮定付きである**(近似で埋めていないことの表明)。 -/
noncomputable def ramificationFiltrationOfCompat
    (hcompat : ∀ K : PAdicLocalField p, ∀ ⦃M⦄, M ∈ openNormalBase K.absGal →
      ∀ ⦃N⦄, N ∈ openNormalBase K.absGal → N ≤ M → ∀ v,
        ((absGalStage K N v : Subgroup K.absGal) : Set K.absGal) * (M : Set K.absGal)
          = ((absGalStage K M v : Subgroup K.absGal) : Set K.absGal)) :
    RamificationFiltration p :=
  ramificationFiltrationOfStages (fun K => absGalStageFiltration K (hcompat K))

end ABC3.Found.PGC
