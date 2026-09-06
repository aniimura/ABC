import ABC3.Found.PGC.UnramifiedGalCharCount
import Mathlib.Topology.Algebra.Category.ProfiniteGrp.Completion
import Mathlib.FieldTheory.Galois.Profinite

/-!
# `Gal(K^ur/K) ≃ₜ* Ẑ`(`sorry` 無し)

経路 Λ の節点 Λ5。Λ4(`Found/PGC/AbelianDecomposition.lean`)が

  `Gal(K_π · K^ur / K) ≃* 𝒪_K^× × Gal(K^ur/K)`

の形で止めていた右因子を、**`Ẑ` と位相群として同定する**。

## `Ẑ` はどこにあるか(2026-09-06 の在庫調査)

* **`Ẑ` は mathlib に在る**——`ProfiniteGrp.ProfiniteCompletion.completion`
  (`Mathlib/Topology/Algebra/Category/ProfiniteGrp/Completion.lean`、Adam Topaz)。
  本ファイルの `ZHat` は `completion (GrpCat.of (Multiplicative ℤ))` である。
  普遍性は `lift` / `lift_eta` / `lift_unique` として在る。
* **`Gal(𝔽̄_q/𝔽_q) ≅ Ẑ` は mathlib に無い**(`.cache/mathlib-index.txt` を
  `ProfiniteCompletion` / `absoluteGaloisGroup` / `GaloisField` / `topologicalGenerator`
  で引いた結果、有限体の絶対 Galois 群を `Ẑ` と同定する宣言は 1 つも無い)。
  在るのは**有限次**の側だけ——`FiniteField.frobeniusAlgEquivOfAlgebraic` と
  `FiniteField.orderOf_frobeniusAlgEquivOfAlgebraic`(Frobenius の位数 = `[L:K]`)。
  したがって「剰余体へ降りて mathlib の結果を引く」道は**無い**。
* `Gal(K^ur/K)` が副有限であること(中間目標 1)は在庫で済む——
  `isGalois_unramifiedClosure`(`UnramifiedExtension.lean`)から
  `InfiniteGalois` の `CompactSpace Gal(K/k)` と `krullTopology` の
  `TotallySeparatedSpace` が付いて `ProfiniteGrp.of` に載る。

## 筋(剰余体を経由しない)

`G := Gal(K^ur/K)`。`UnramifiedGalCharCount.lean` が用意している
`exists_unramified_surjective_zmod`(`G ↠ ℤ/N` で核が次数 `N` の不分岐中間体の
固定部分群)と `exists_unramified_fixingSubgroup_le`(開部分群はどれかの
`P_N` を含む)を、**選択して名前を付ける**ところから始める。

1. **段** `unramLevel K N`(`K^ur` の中の次数 `N` の不分岐中間体)と
   `unramLevelHom K N : G ↠ ℤ/N`(核 = `P_N := (unramLevel K N).fixingSubgroup`)。
   `adjoin_eq_of_isUnramified`(同じ次数の不分岐拡大は一意)により
   `unramLevel` は生成元の取り方に依らず、`M ∣ N ⟹ P_N ≤ P_M` が従う。
2. **整合的 Frobenius**(`exists_coherentFrobenius`)。
   `unramLevelGeneratorSet K N := {g | ∀ h, ∃ k, h⁻¹ gᵏ ∈ P_N}` は
   `unramLevelHom K N` の逆像なので**開かつ閉**(核が開な準同型の逆像は
   clopen——`isClosed_preimage_of_isOpen_ker`)、空でなく、`N ∣ NM` から
   有向に減少する。`G` はコンパクトなので
   `IsCompact.nonempty_iInter_of_directed_nonempty_isCompact_isClosed` で
   **全段を同時に生成する `σ`** が取れる。
   ★これが「Frobenius が位相的生成元」の中身であり、剰余体の言葉を一度も
   使わずに済ませている点が本ファイルの要である。
3. **同型**。`φ : Multiplicative ℤ →* G`(`1 ↦ σ`)に普遍性を当てて
   `ψ := lift φ : Ẑ ⟶ G` を作り、全単射を示す。
   * 全射:`range ψ` はコンパクト(⟹ 閉)で `range φ` を含み、
     `range φ` は稠密(`dense_zpowers_frobenius`——開近傍は `P_N` を含み、
     `σ` は `G/P_N` を生成する)。
   * 単射:`lift φ ≫ proj U = π_{φ⁻¹U} ≫ quotientMap φ U`(`lift_comp_proj`、
     `lift_unique` で 1 行)。`quotientMap φ U` は常に単射なので
     `ψ x = 1 ⟹ x_{φ⁻¹U} = 1`。`H` を指数 `N` の有限指数正規部分群とすると
     `φ⁻¹(P_N) ≤ H`(`σ` の `G/P_N` での位数がちょうど `N` だから
     `σᵏ ∈ P_N ⟺ N ∣ k`、そして `Subgroup.pow_index_mem`)なので
     極限の整合性から `x_H = 1`。

## ★設計上の注意(守ったこと)

* **`inertia` を経由していない**。不分岐側は `unramifiedClosure K` という
  体として直接扱っており、`Interface` の `SubgroupCorrespondence` /
  `ResidueCardinality` は本ファイルの主張にも証明にも現れない
  (Corollary 1.3 との循環を避けるため)。
  ★ただしファイル依存としては `UnramifiedGalCharCount` の import 閉包に
  既存構造が入っている(`AbelianDecomposition.lean` と同じ状況)。
* **`Abelianization` を使っていない**(副有限群では有限指数部分群が開とは
  限らないため)。商群は `ProfiniteGrp` の極限記述の中でしか作らない。
* **結論に自由なパラメータを出していない**——
  `unramifiedClosureGalEquivZHat K` の型は `K` にしか依存しない。
  `σ` は `exists_coherentFrobenius` の `Classical.choose` に閉じ込めてある。

## 逸脱(記録)

* 古典的な証明は `Gal(K^ur/K) ≅ Gal(𝔽̄_q/𝔽_q)`(不分岐の定義)を経由して
  有限体の Frobenius を使う。本ファイルは**剰余体へ降りない**——
  「次数 `N` の不分岐拡大が各 `N` にちょうど 1 つあり、その Galois 群が
  巡回群 `ℤ/N` である」(`UnramifiedExtension.lean` が既に持っている)
  という事実だけから、生成元をコンパクト性で束ねる。
  したがって「Frobenius」という名前は本ファイルでは
  **`exists_coherentFrobenius` が返す元**を指し、剰余体上の `x ↦ x^q` との
  一致は主張していない(それは `UnramifiedFrobenius.lean` の担当)。
* `Ẑ` は `Multiplicative ℤ` の副有限完備化として与える(mathlib の
  `ProfiniteCompletion` が乗法群の圏 `GrpCat` の上にあるため)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open CategoryTheory ProfiniteGrp ProfiniteGrp.ProfiniteCompletion

variable {p : ℕ} [Fact p.Prime]

/-- `Gal(K^ur/K)`。 -/
abbrev unramGal (K : PAdicLocalField p) : Type :=
  ↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K)

instance instIsGaloisUnramifiedClosure (K : PAdicLocalField p) :
    IsGalois K.carrier ↥(unramifiedClosure K) := isGalois_unramifiedClosure K

/-! ## 1. 不分岐の塔に名前を付ける

`UnramifiedGalCharCount.lean::exists_unramified_surjective_zmod` は
「次数 `N` の不分岐中間体と `G ↠ ℤ/N` が在る」ことしか言わない。
`Ẑ` を作るには**段どうしの整合性**が要るので、選択して名前を付け、
`adjoin_eq_of_isUnramified`(一意性)で整合性を回復する。 -/

/-- `K^ur` の中の次数 `N` の不分岐中間体の生成元。`N = 0` では `0`(意味を持たない)。 -/
noncomputable def unramLevelGen (K : PAdicLocalField p) (N : ℕ) : ↥(unramifiedClosure K) :=
  if h : N = 0 then 0 else (exists_unramified_surjective_zmod K h).choose

/-- `K^ur` の中の次数 `N` の不分岐中間体。 -/
noncomputable def unramLevel (K : PAdicLocalField p) (N : ℕ) :
    IntermediateField K.carrier ↥(unramifiedClosure K) :=
  IntermediateField.adjoin K.carrier {unramLevelGen K N}

/-- `Gal(K^ur/K) ↠ ℤ/N`(核は `unramLevel K N` の固定部分群)。 -/
noncomputable def unramLevelHom (K : PAdicLocalField p) (N : ℕ) :
    unramGal K →* Multiplicative (ZMod N) :=
  if h : N = 0 then 1 else ((exists_unramified_surjective_zmod K h).choose_spec).choose

theorem unramLevelGen_of_ne (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0) :
    unramLevelGen K N = (exists_unramified_surjective_zmod K hN).choose := dif_neg hN

theorem unramLevelHom_of_ne (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0) :
    unramLevelHom K N = ((exists_unramified_surjective_zmod K hN).choose_spec).choose := dif_neg hN

/-- `unramLevel` / `unramLevelHom` の仕様(`N ≠ 0`)。

退化の自己検査:`N = 0` では `ZMod 0 = ℤ` で全射になりえない
(`Gal(K^ur/K)` の有限商は `ℤ` を持たない)ので `hN` は落とせない。 -/
theorem unramLevel_spec (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0) :
    IsUnramifiedAdjoin K ((unramLevelGen K N : K.closure)) ∧
    Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier
        ({((unramLevelGen K N : ↥(unramifiedClosure K)) : K.closure)} : Set K.closure)) = N ∧
    Function.Surjective (unramLevelHom K N) ∧
    (unramLevelHom K N).ker = (unramLevel K N).fixingSubgroup := by
  rw [unramLevel, unramLevelGen_of_ne K hN, unramLevelHom_of_ne K hN]
  exact ((exists_unramified_surjective_zmod K hN).choose_spec).choose_spec

instance finiteDimensional_unramLevel (K : PAdicLocalField p) (N : ℕ) :
    FiniteDimensional K.carrier (unramLevel K N) :=
  IntermediateField.finiteDimensional_adjoin (fun z _ => Algebra.IsIntegral.isIntegral z)

theorem finrank_unramLevel (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0) :
    Module.finrank K.carrier (unramLevel K N) = N :=
  ((adjoinUnramifiedAlgEquiv K (unramLevelGen K N)).toLinearEquiv.finrank_eq).trans
    (unramLevel_spec K hN).2.1

/-- **塔**:`M ∣ N ⟹ `unramLevel K M ≤ unramLevel K N`。 -/
theorem unramLevel_le_unramLevel_of_dvd (K : PAdicLocalField p) {M N : ℕ} (hM : M ≠ 0)
    (hN : N ≠ 0) (h : M ∣ N) : unramLevel K M ≤ unramLevel K N :=
  adjoin_le_adjoin_of_val_le K _ _
    (adjoin_le_of_dvd K _ _ (unramLevel_spec K hM).1 (unramLevel_spec K hN).1
      (by rw [(unramLevel_spec K hM).2.1, (unramLevel_spec K hN).2.1]; exact h))

/-- 塔の包含は固定部分群の逆包含。 -/
theorem fixingSubgroup_unramLevel_le (K : PAdicLocalField p) {M N : ℕ} (hM : M ≠ 0) (hN : N ≠ 0)
    (h : M ∣ N) : (unramLevel K N).fixingSubgroup ≤ (unramLevel K M).fixingSubgroup :=
  IntermediateField.fixingSubgroup_le (unramLevel_le_unramLevel_of_dvd K hM hN h)

/-- **一意性**:次数 `N` の不分岐中間体は `unramLevel K N` に限る。 -/
theorem unramLevel_eq_adjoin (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0)
    (y : ↥(unramifiedClosure K)) (hy : IsUnramifiedAdjoin K (y : K.closure))
    (hrank : Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({(y : K.closure)} : Set K.closure)) = N) :
    unramLevel K N = IntermediateField.adjoin K.carrier ({y} : Set ↥(unramifiedClosure K)) := by
  have heq : IntermediateField.adjoin K.carrier
      ({((unramLevelGen K N : ↥(unramifiedClosure K)) : K.closure)} : Set K.closure)
      = IntermediateField.adjoin K.carrier ({(y : K.closure)} : Set K.closure) :=
    adjoin_eq_of_isUnramified K _ _ (unramLevel_spec K hN).1 hy
      (by rw [(unramLevel_spec K hN).2.1, hrank])
  exact le_antisymm (adjoin_le_adjoin_of_val_le K _ _ heq.le)
    (adjoin_le_adjoin_of_val_le K _ _ heq.ge)

theorem isOpen_fixingSubgroup_unramLevel (K : PAdicLocalField p) (N : ℕ) :
    IsOpen (((unramLevel K N).fixingSubgroup : Subgroup (unramGal K)) : Set (unramGal K)) :=
  IntermediateField.fixingSubgroup_isOpen _

/-- **開部分群はどれかの `P_N` を含む**。

退化の自己検査:`hU`(開)を落とすと**偽**——閉でない部分群や、有限指数でも
開でない部分群には基本近傍が入らない。 -/
theorem exists_unramLevel_fixingSubgroup_le (K : PAdicLocalField p) {U : Subgroup (unramGal K)}
    (hU : IsOpen (U : Set (unramGal K))) :
    ∃ N : ℕ, N ≠ 0 ∧ (unramLevel K N).fixingSubgroup ≤ U := by
  obtain ⟨y, hyu, hyle⟩ := exists_unramified_fixingSubgroup_le K hU
  refine ⟨Module.finrank K.carrier
    (IntermediateField.adjoin K.carrier ({(y : K.closure)} : Set K.closure)), ?_, ?_⟩
  · exact Module.finrank_pos.ne'
  · rw [unramLevel_eq_adjoin K Module.finrank_pos.ne' y hyu rfl]
    exact hyle

theorem index_fixingSubgroup_unramLevel (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0) :
    ((unramLevel K N).fixingSubgroup : Subgroup (unramGal K)).index = N :=
  (IntermediateField.finrank_eq_fixingSubgroup_index (unramLevel K N)).symm.trans
    (finrank_unramLevel K hN)

theorem isGalois_unramLevel (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0) :
    IsGalois K.carrier (unramLevel K N) := by
  haveI : Normal K.carrier
      (IntermediateField.adjoin K.carrier
        ({((unramLevelGen K N : ↥(unramifiedClosure K)) : K.closure)} : Set K.closure)) :=
    normal_of_isUnramifiedAdjoin K _ (unramLevel_spec K hN).1
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  haveI : Algebra.IsSeparable K.carrier
      (IntermediateField.adjoin K.carrier
        ({((unramLevelGen K N : ↥(unramifiedClosure K)) : K.closure)} : Set K.closure)) :=
    IntermediateField.isSeparable_tower_bot K.carrier _
  haveI : IsGalois K.carrier
      (IntermediateField.adjoin K.carrier
        ({((unramLevelGen K N : ↥(unramifiedClosure K)) : K.closure)} : Set K.closure)) := ⟨⟩
  exact IsGalois.of_algEquiv (adjoinUnramifiedAlgEquiv K (unramLevelGen K N)).symm

theorem normal_fixingSubgroup_unramLevel (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0) :
    ((unramLevel K N).fixingSubgroup : Subgroup (unramGal K)).Normal :=
  (InfiniteGalois.normal_iff_isGalois (unramLevel K N)).mpr (isGalois_unramLevel K hN)

/-! ## 2. 整合的 Frobenius(コンパクト性)

各段の生成元は在るが、**段をまたいで整合する**生成元は選択公理だけでは
取れない。生成元の集合が clopen であることと `G` のコンパクト性を使う。 -/

/-- 核が開な準同型による逆像は開。 -/
theorem isOpen_preimage_of_isOpen_ker {G A : Type*} [Group G] [TopologicalSpace G]
    [IsTopologicalGroup G] [Group A] (f : G →* A)
    (hf : IsOpen ((f.ker : Subgroup G) : Set G)) (S : Set A) : IsOpen (f ⁻¹' S) := by
  rw [isOpen_iff_forall_mem_open]
  intro x hx
  refine ⟨(fun g => x * g) '' ((f.ker : Subgroup G) : Set G), ?_,
    (isOpenMap_mul_left x) _ hf, ⟨1, one_mem _, mul_one x⟩⟩
  rintro _ ⟨k, hk, rfl⟩
  simpa only [Set.mem_preimage, map_mul, MonoidHom.mem_ker.mp hk, mul_one] using hx

/-- 核が開な準同型による逆像は閉。 -/
theorem isClosed_preimage_of_isOpen_ker {G A : Type*} [Group G] [TopologicalSpace G]
    [IsTopologicalGroup G] [Group A] (f : G →* A)
    (hf : IsOpen ((f.ker : Subgroup G) : Set G)) (S : Set A) : IsClosed (f ⁻¹' S) := by
  rw [← isOpen_compl_iff, ← Set.preimage_compl]
  exact isOpen_preimage_of_isOpen_ker f hf _

/-- 像が `Gal(K^ur/K)/P_N` を生成するような元の全体。 -/
def unramLevelGeneratorSet (K : PAdicLocalField p) (N : ℕ) : Set (unramGal K) :=
  {g | ∀ h : unramGal K, ∃ k : ℕ, h⁻¹ * g ^ k ∈ (unramLevel K N).fixingSubgroup}

theorem unramLevelGeneratorSet_eq_preimage (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0) :
    unramLevelGeneratorSet K N
      = unramLevelHom K N ⁻¹' {a : Multiplicative (ZMod N) | ∀ b, ∃ k : ℕ, a ^ k = b} := by
  obtain ⟨-, -, hsurj, hker⟩ := unramLevel_spec K hN
  ext g
  simp only [unramLevelGeneratorSet, Set.mem_setOf_eq, Set.mem_preimage]
  constructor
  · intro hg b
    obtain ⟨h, rfl⟩ := hsurj b
    obtain ⟨k, hk⟩ := hg h
    rw [← hker, MonoidHom.mem_ker, map_mul, map_inv, map_pow] at hk
    exact ⟨k, (inv_mul_eq_one.mp hk).symm⟩
  · intro hg h
    obtain ⟨k, hk⟩ := hg (unramLevelHom K N h)
    refine ⟨k, ?_⟩
    rw [← hker, MonoidHom.mem_ker, map_mul, map_inv, map_pow, hk, inv_mul_cancel]

theorem isClosed_unramLevelGeneratorSet (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0) :
    IsClosed (unramLevelGeneratorSet K N) := by
  rw [unramLevelGeneratorSet_eq_preimage K hN]
  exact isClosed_preimage_of_isOpen_ker _
    (by rw [(unramLevel_spec K hN).2.2.2]; exact isOpen_fixingSubgroup_unramLevel K N) _

theorem unramLevelGeneratorSet_nonempty (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0) :
    (unramLevelGeneratorSet K N).Nonempty := by
  haveI : NeZero N := ⟨hN⟩
  obtain ⟨-, -, hsurj, -⟩ := unramLevel_spec K hN
  obtain ⟨g, hg⟩ := hsurj (Multiplicative.ofAdd (1 : ZMod N))
  refine ⟨g, ?_⟩
  rw [unramLevelGeneratorSet_eq_preimage K hN]
  intro b
  exact ⟨(Multiplicative.toAdd b).val, by rw [hg, ofAdd_one_pow_val]; rfl⟩

theorem unramLevelGeneratorSet_mono (K : PAdicLocalField p) {M N : ℕ} (hM : M ≠ 0) (hN : N ≠ 0)
    (h : M ∣ N) : unramLevelGeneratorSet K N ⊆ unramLevelGeneratorSet K M := by
  intro g hg h'
  obtain ⟨k, hk⟩ := hg h'
  exact ⟨k, fixingSubgroup_unramLevel_le K hM hN h hk⟩

/-- **★★★★★★★★★★すべての段を同時に生成する元(Frobenius)が存在する**。

各 `N` に対する生成元の集合は空でない閉集合で、`N ∣ NM` により有向に減少する。
`Gal(K^ur/K)` はコンパクトなので共通部分は空でない。

★ここが「Frobenius が位相的生成元」の中身である。剰余体の言葉は使っていない。 -/
theorem exists_coherentFrobenius (K : PAdicLocalField p) :
    ∃ σ : unramGal K, ∀ N : ℕ, N ≠ 0 → σ ∈ unramLevelGeneratorSet K N := by
  haveI : Nonempty {N : ℕ // N ≠ 0} := ⟨⟨1, one_ne_zero⟩⟩
  obtain ⟨σ, hσ⟩ := IsCompact.nonempty_iInter_of_directed_nonempty_isCompact_isClosed
    (fun i : {N : ℕ // N ≠ 0} => unramLevelGeneratorSet K i.1)
    (fun i j => ⟨⟨i.1 * j.1, mul_ne_zero i.2 j.2⟩,
      unramLevelGeneratorSet_mono K i.2 (mul_ne_zero i.2 j.2) ⟨j.1, rfl⟩,
      unramLevelGeneratorSet_mono K j.2 (mul_ne_zero i.2 j.2) ⟨i.1, mul_comm i.1 j.1⟩⟩)
    (fun i => unramLevelGeneratorSet_nonempty K i.2)
    (fun i => (isClosed_unramLevelGeneratorSet K i.2).isCompact)
    (fun i => isClosed_unramLevelGeneratorSet K i.2)
  exact ⟨σ, fun N hN => Set.mem_iInter.mp hσ ⟨N, hN⟩⟩

/-- 生成元の `G/P_N` での位数はちょうど `N`。 -/
theorem orderOf_unramLevelHom_of_mem (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0)
    {σ : unramGal K} (hσ : σ ∈ unramLevelGeneratorSet K N) :
    orderOf (unramLevelHom K N σ) = N := by
  haveI : NeZero N := ⟨hN⟩
  rw [unramLevelGeneratorSet_eq_preimage K hN] at hσ
  have htop : Subgroup.zpowers (unramLevelHom K N σ) = ⊤ := by
    rw [Subgroup.eq_top_iff']
    intro b
    obtain ⟨k, hk⟩ := hσ b
    exact ⟨(k : ℤ), by simpa [zpow_natCast] using hk⟩
  rw [← Nat.card_zpowers, htop]
  simp [Nat.card_eq_fintype_card]

/-- **`σᵏ ∈ P_N ⟺ N ∣ k`**。 -/
theorem zpow_mem_fixingSubgroup_unramLevel_iff (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0)
    {σ : unramGal K} (hσ : σ ∈ unramLevelGeneratorSet K N) (k : ℤ) :
    σ ^ k ∈ (unramLevel K N).fixingSubgroup ↔ (N : ℤ) ∣ k := by
  rw [← (unramLevel_spec K hN).2.2.2, MonoidHom.mem_ker, map_zpow,
    ← orderOf_dvd_iff_zpow_eq_one, orderOf_unramLevelHom_of_mem K hN hσ]

/-! ## 3. `Ẑ` との同型 -/

/-- `Gal(K^ur/K)` を副有限群として。 -/
noncomputable abbrev unramGalProfinite (K : PAdicLocalField p) : ProfiniteGrp.{0} :=
  ProfiniteGrp.of (unramGal K)

/-- **`Ẑ`** —— `ℤ` の副有限完備化(mathlib の `ProfiniteCompletion`)。 -/
noncomputable abbrev ZHat : ProfiniteGrp.{0} :=
  ProfiniteGrp.ProfiniteCompletion.completion (GrpCat.of (Multiplicative ℤ))

/-- `1 ↦ σ` で定まる `ℤ → Gal(K^ur/K)`。 -/
noncomputable def frobGrpHom (K : PAdicLocalField p) (σ : unramGal K) :
    GrpCat.of (Multiplicative ℤ) ⟶ GrpCat.of (unramGalProfinite K) :=
  GrpCat.ofHom (zpowersHom (unramGal K) σ)

/-- 副有限完備化の普遍性で持ち上げた `Ẑ → Gal(K^ur/K)`。 -/
noncomputable def frobLift (K : PAdicLocalField p) (σ : unramGal K) :
    ZHat ⟶ unramGalProfinite K :=
  ProfiniteGrp.ProfiniteCompletion.lift (frobGrpHom K σ)

/-- `P_N` を開正規部分群として。 -/
noncomputable def unramLevelOpenNormal (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0) :
    OpenNormalSubgroup (unramGalProfinite K) where
  toOpenSubgroup := ⟨(unramLevel K N).fixingSubgroup, isOpen_fixingSubgroup_unramLevel K N⟩
  isNormal' := normal_fixingSubgroup_unramLevel K hN

/-- 普遍射を開正規部分群の商へ落とすと、`Ẑ` の該当成分に `quotientMap` を当てたもの。

`lift_unique` で 1 行(両辺は `eta` を前置すると `y ↦ mk (f y)` で一致する)。 -/
theorem lift_comp_proj {G : GrpCat.{0}} {P : ProfiniteGrp.{0}} (f : G ⟶ GrpCat.of P)
    (U : OpenNormalSubgroup P) :
    lift f ≫ ProfiniteGrp.proj U
      = (ProfiniteGrp.limitCone (diagram G)).π.app (preimage f U)
        ≫ ProfiniteGrp.ofFiniteGrpHom (quotientMap f U) := by
  apply lift_unique
  rw [Functor.map_comp, ← Category.assoc, lift_eta]
  ext y
  rfl

/-- `quotientMap f U : G/f⁻¹U → P/U` は常に単射。 -/
theorem eq_one_of_quotientMap_eq_one {G : GrpCat.{0}} {P : ProfiniteGrp.{0}}
    (f : G ⟶ GrpCat.of P) (U : OpenNormalSubgroup P)
    (z : G ⧸ (preimage f U).toSubgroup) (hz : quotientMap f U z = 1) : z = 1 := by
  induction z using QuotientGroup.induction_on with
  | H a =>
    rw [QuotientGroup.eq_one_iff]
    exact (QuotientGroup.eq_one_iff (f.hom a)).mp hz

/-- **`⟨σ⟩` は `Gal(K^ur/K)` で稠密**。

`1` の開近傍は Krull 位相の基本近傍 `E.fixingSubgroup` を含み、それは
どれかの `P_N` を含む(`exists_unramLevel_fixingSubgroup_le`)。 -/
theorem dense_zpowers_frobenius (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ N : ℕ, N ≠ 0 → σ ∈ unramLevelGeneratorSet K N) :
    Dense (Set.range (fun n : ℤ => σ ^ n)) := by
  rw [dense_iff_inter_open]
  rintro U hU ⟨x, hx⟩
  have hcont : Continuous (fun y : unramGal K => x * y) := continuous_const.mul continuous_id
  have hnhd : {y : unramGal K | x * y ∈ U} ∈ nhds (1 : unramGal K) :=
    (hU.preimage hcont).mem_nhds (by simpa using hx)
  obtain ⟨E, hEfin, hEsub⟩ :=
    (krullTopology_mem_nhds_one_iff K.carrier ↥(unramifiedClosure K) _).mp hnhd
  haveI := hEfin
  obtain ⟨N, hN, hNle⟩ := exists_unramLevel_fixingSubgroup_le K
    (IntermediateField.fixingSubgroup_isOpen E)
  obtain ⟨k, hk⟩ := hσ N hN x
  refine ⟨σ ^ (k : ℤ), ?_, ⟨(k : ℤ), rfl⟩⟩
  have hmem : x * (x⁻¹ * σ ^ k) ∈ U := hEsub (hNle hk)
  rwa [mul_inv_cancel_left, ← zpow_natCast] at hmem

/-- 全射:像はコンパクト(⟹ 閉)で、稠密な `⟨σ⟩` を含む。 -/
theorem frobLift_surjective (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ N : ℕ, N ≠ 0 → σ ∈ unramLevelGeneratorSet K N) :
    Function.Surjective (frobLift K σ) := by
  have hcont : Continuous (frobLift K σ) := (frobLift K σ).hom.continuous_toFun
  have hclosed : IsClosed (Set.range (frobLift K σ)) := (isCompact_range hcont).isClosed
  have hsub : Set.range (fun n : ℤ => σ ^ n) ⊆ Set.range (frobLift K σ) := by
    rintro _ ⟨n, rfl⟩
    exact ⟨ProfiniteCompletion.etaFn _ (Multiplicative.ofAdd n),
      ConcreteCategory.congr_hom (ProfiniteCompletion.lift_eta (frobGrpHom K σ))
        (Multiplicative.ofAdd n)⟩
  rw [← Set.range_eq_univ]
  refine Set.eq_univ_of_univ_subset ?_
  rw [← (dense_zpowers_frobenius K hσ).closure_eq, ← hclosed.closure_eq]
  exact closure_mono hsub

/-- 単射:`ψ x = 1` なら、指数 `N` の有限指数正規部分群 `H` に対し
`φ⁻¹(P_N) ≤ H` なので極限の整合性から `x_H = 1`。 -/
theorem frobLift_injective (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ N : ℕ, N ≠ 0 → σ ∈ unramLevelGeneratorSet K N) :
    Function.Injective (frobLift K σ) := by
  have key : ∀ x : ZHat, frobLift K σ x = 1 → x = 1 := by
    intro x hx
    refine ProfiniteGrp.limit_ext _ x 1 (fun H => ?_)
    rw [ProfiniteGrp.limit_one_val]
    haveI := H.isNormal'
    have hN : H.toSubgroup.index ≠ 0 := Subgroup.finiteIndex_iff.mp H.isFiniteIndex'
    have hle : preimage (frobGrpHom K σ) (unramLevelOpenNormal K hN) ≤ H := by
      intro a ha
      obtain ⟨n, rfl⟩ : ∃ n : ℤ, (Multiplicative.ofAdd n : Multiplicative ℤ) = a :=
        ⟨Multiplicative.toAdd (show Multiplicative ℤ from a), rfl⟩
      have ha' : σ ^ n ∈ (unramLevel K H.toSubgroup.index).fixingSubgroup := ha
      obtain ⟨m, hm⟩ := (zpow_mem_fixingSubgroup_unramLevel_iff K hN (hσ _ hN) n).mp ha'
      have hpow := H.toSubgroup.pow_index_mem (Multiplicative.ofAdd m)
      have heq : (Multiplicative.ofAdd m) ^ H.toSubgroup.index
          = (Multiplicative.ofAdd n : Multiplicative ℤ) := by
        rw [← ofAdd_nsmul, nsmul_eq_mul, ← hm]
      exact heq ▸ hpow
    have hx' : (lift (frobGrpHom K σ)) x = 1 := hx
    have hLHS : (ProfiniteGrp.proj (unramLevelOpenNormal K hN))
        ((lift (frobGrpHom K σ)) x) = 1 := by
      rw [hx']; exact map_one _
    have h2 := ConcreteCategory.congr_hom
      (lift_comp_proj (frobGrpHom K σ) (unramLevelOpenNormal K hN)) x
    simp only [ProfiniteGrp.comp_apply] at h2
    have h1 : x.val (preimage (frobGrpHom K σ) (unramLevelOpenNormal K hN)) = 1 :=
      eq_one_of_quotientMap_eq_one _ _ _ (h2.symm.trans hLHS)
    have hcompat := x.2 hle.hom
    rw [← hcompat, h1]
    exact map_one _
  intro a b hab
  have hone : (frobLift K σ) (a * b⁻¹) = 1 := by
    rw [map_mul, map_inv, hab, mul_inv_cancel]
  have := key _ hone
  rwa [mul_inv_eq_one] at this

/-- `Ẑ ≃ₜ* Gal(K^ur/K)`(生成元 `σ` を明示した版)。 -/
noncomputable def zhatContinuousMulEquivUnramGal (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ N : ℕ, N ≠ 0 → σ ∈ unramLevelGeneratorSet K N) : ZHat ≃ₜ* unramGal K :=
  { Continuous.homeoOfEquivCompactToT2
      (f := Equiv.ofBijective _ ⟨frobLift_injective K hσ, frobLift_surjective K hσ⟩)
      (frobLift K σ).hom.continuous_toFun with
    map_mul' := (frobLift K σ).hom.map_mul' }

/-- **★★★★★★★★★★★★★★★★★★★★(Λ5)`Gal(K^ur/K) ≃ₜ* Ẑ`**。

経路 Λ の節点 Λ5。`Ẑ = ProfiniteGrp.ProfiniteCompletion.completion (GrpCat.of (Multiplicative ℤ))`。

**証明**:整合的 Frobenius `σ`(`exists_coherentFrobenius`)を取り、
副有限完備化の普遍性で `ψ : Ẑ → Gal(K^ur/K)` を作る。`ψ` は
稠密像 + コンパクト性で全射、極限の成分計算で単射。

退化の自己検査。

* `K^ur` を `K.closure` に替えると**偽**——`Γ_K` は非可換
  (`QpNonAbelian.lean`)なので `Ẑ` とは同型になりえない。
* `K^ur` を有限次不分岐拡大 `K_N` に替えると**偽**(左辺が有限群 `ℤ/N`)。
* 位相を落として抽象群の同型として読むと、主張は依然として正しいが
  本ファイルの証明は使えない(全射性がコンパクト性に依存している)。

★結論の型は `K` にしか依存しない——`σ` は `Classical.choose` に閉じ込めてある。 -/
noncomputable def unramifiedClosureGalEquivZHat (K : PAdicLocalField p) :
    (↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K)) ≃ₜ* ZHat :=
  (zhatContinuousMulEquivUnramGal K (exists_coherentFrobenius K).choose_spec).symm

/-- `Gal(K^ur/K) ≃ₜ* Ẑ`(命題の形)。 -/
theorem nonempty_unramifiedClosureGalEquivZHat (K : PAdicLocalField p) :
    Nonempty ((↥(unramifiedClosure K) ≃ₐ[K.carrier] ↥(unramifiedClosure K)) ≃ₜ* ZHat) :=
  ⟨unramifiedClosureGalEquivZHat K⟩

end ABC3.Found.PGC
