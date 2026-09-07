import ABC3.Found.PGC.LocalClassFieldTheory
import ABC3.Found.PGC.LubinTateClosureTopology
import ABC3.Found.PGC.UnitsPowP

/-!
# Artin 写像 `Art_π` —— Yoshida 2008 §4.2 / Definition 4.10

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Section 4.2(物理 p.7)と
Definition 4.10(物理 p.9)。構造化済み原文は
`ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/section-4.html`
の `#setup-4-2` と `#def-4-10`。

原文 (Yoshida08 p.9):
> Definition 4.10. For any f ∈ O[scr]_L[X] with L/K finite, set K^m := K^urL^m_f. Then K^m/K is finitely ramified, and Galois by Proposition 4.7(i). By Lemma 2.2, the completion of K^m is KL^m_f = K^m and K^m = K^m ∩ K^sep, thus independent of f. Setting K^LT := _m≥1 K^m = K^LT ∩K^sep, we have W(K^LT/K) ∼ = W(K^LT/K) by the remark after Definition 2.5. We call a finite extension of K a Lubin-Tate extension if it is contained in K^LT. We call the inverse of ρ the Artin map of K and write Art_K : K^× ∼ =−→ W(K^LT/K). We have v ◦ Art_K = v.

## 何を作ったか

`Art_π : K^× →* Gal(K^ab/K)` を **単項準同型として実際に構成する**。骨格は 2 つの直積分解:

| 側 | 分解 | 在庫 |
|---|---|---|
| 源 `K^×` | `K^× ≅ ℤ × 𝒪_K^×`（π で分裂） | `UnitsSplit.lean::unitsSplitEquiv` |
| 標的 `Gal(K^ab/K)` | `≅ 𝒪_K^× × Gal(K^ur/K)` | `AbelianDecomposition.lean::lubinTateUnramifiedGalEquivProd` + **LKW** |

★標的側の分解が `K^ab` について言えるのは、**Theorem 6.15（局所 Kronecker-Weber）が
`K^ab = K_π · K^ur` を与えるから**である（`LocalClassFieldTheory.lean`）。
本ファイルはそれを `abelianClosure_eq_lubinTateClosure_sup_unramifiedClosure` の形で
**実際に消費している**（§6 `artinMapAbelian`）。LKW が無ければ `Gal(K^ab/K)` の分解は
書けず、Artin 写像の標的は `Gal(K_π·K^ur/K)` に留まる。

`Art_π` の値は原典 §4.2 の正規化にそのまま従う（下記「符号」）:

* `u ∈ 𝒪_K^×` 上: `Art_π(u)` は `K_π` 上で相互律同型 `lubinTateClosureGalEquivUnits` の
  **逆**（＝Lubin-Tate 加群への `[u]` 倍作用）、`K^ur` 上では恒等。
* `π` 上: `Art_π(π)` は `K_π` 上で恒等、`K^ur` 上で**幾何 Frobenius**
  `Frob_K = ϕ^{-1}`（`geometricFrobenius`）。

★**両方を落としていない**（片方だけでは写像が決まらない）:
`artinMap_unitsToCarrier`（`𝒪^×` 上）と `artinMap_uniformizer`（`π` 上）の 2 本が
それぞれ値を決めており、`K^× = π^ℤ · 𝒪^×` なのでこの 2 本で `Art_π` は一意に決まる。

## ★符号（原典の正規化は幾何 Frobenius である）

原典 §2.2 は Weil 群を `W(E′/K) := {σ | σ|_E ∈ Frob_K^ℤ}` と定義し、`Frob_K := ϕ^{-1}`
（幾何 Frobenius、`ϕ` は算術 Frobenius）と置く。Proposition 4.7(ii) の
`ρ_{f,m}` は「`K̂` 上 `ϕ^j`、`µ_{f,m}` 上 `α ↦ [xπ_j](α)`」なる `σ` を
`x mod (1+𝔭^m)`（ただし `v(x) = −j`）へ送るから、`v(x) = 1` の `x`（＝素元）に
対応する `σ` は `K̂` 上 `ϕ^{-1} = Frob_K` である。ゆえに Definition 4.10 の
`v ∘ Art_K = v` は「`Art_K(x)|_{K^ur} = Frob_K^{v(x)}`」の意味になる。
★本ファイルはこの正規化を採る（`geometricFrobenius`）。
`artinMap_unramified_component` が `v ∘ Art = v` にあたる。

## ★★射程 —— 何を主張していないか

* ★★**`Art_π` は `π`（と Lubin-Tate 級数 `f`）に依存する。**
  ★**「Art は π に依らない」とは主張していない。** 非依存性は原典 Corollary 4.9
  （`ρ_{f,m}` の `f` 非依存）＋ Dwork の定理を要する**別のノード**であり、
  本ファイルの射程外である（同時走行の `LubinTateUniformizerIndependence` が担当）。
  ★依存は型に現れている: `artinMap` の引数に `hπmax`・`f`・`hf0`・`hf1`・`hf` がある。
* ★★**全射だとは主張していない。** `Gal(K^ur/K) ≅ Ẑ` に対して像の第 2 成分は
  `Frob_K^ℤ` にすぎない。`range_artinMap` が像を
  **ちょうど Weil 群 `W(K^ab/K) = {σ | σ|_{K^ur} ∈ Frob^ℤ}`** と同定しており、
  これが原典 Definition 4.10 の `Art_K : K^× ≅ W(K^LT/K)` に対応する。
* ★★**像は `Gal(K^ab/K)` の中で稠密**（`dense_range_artinMap`、★仮定なし）。
  `Φ` が**位相群の同型**であること（`continuous_abelianGalEquivProd` /
  `_symm`）を経由する。★稠密であって**全射ではない**。
* ★本ファイルは Proposition 4.7(ii) の `ρ` そのもの（Weil 群からの同型）は作らない。
  `WeilReciprocityExtension.lean` の抽象核がそれを持っており、ここで作るのは
  **`𝒪^×` 側と `ℤ` 側の値を指定して組み上げた具体的な準同型**である。

## 設計 —— 抽象核と具体層

★§1 の抽象核 8 本には**分岐・付値・Galois・Lubin-Tate の語彙が 1 つも出てこない**
（`splitTransportHom` 系 6 本は純群論、`dense_of_dense_image` は位相だけ）。
§2 の 2 本（`autCongr_equivOfEq_restrictNormalHom` /
`galSupEquivProd_restrictNormalHom`）は一般の体と Galois だけで書ける。

★#59（中間体 2 層をまたぐ `rfl`）は**定型 (b)/(c) で回避した**:
`Gal(K^ab/K)` と `Gal(K_π/K)`・`Gal(K^ur/K)` を直接繋がず、
**すべて `Γ_K = Gal(K̄/K)` からの `restrictNormalHom`（1 層）で書いている**
（`abelianGalEquivProd_restrictNormalHom`）。中間体の中の中間体は 1 度も作らない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-! ## §1 抽象核 —— 純群論と位相だけ

★分岐・付値・Galois・Lubin-Tate の語彙が 1 つも出てこない。 -/

section AbstractCore

variable {Z U S A B : Type*} [CommGroup Z] [CommGroup U] [CommGroup S] [Group A] [Group B]

/-- ★★**抽象核 1** —— 「2 つの直積分解のあいだの準同型」。

源が `Z × U ≃* S`、標的が `A ≃* U × B` と分解していて、`ψ : Z →* B` が与えられたとき、

    S ≃ Z × U → U × B ≃ A,   (z, u) ↦ (u, ψ z)

は準同型 `S →* A` である。★これが Artin 写像の骨格そのもので、
局所体・分岐・Lubin-Tate は 1 つも要らない。 -/
def splitTransportHom (e : Z × U ≃* S) (Φ : A ≃* U × B) (ψ : Z →* B) : S →* A :=
  (Φ.symm.toMonoidHom).comp
    (((MonoidHom.snd Z U).prod (ψ.comp (MonoidHom.fst Z U))).comp e.symm.toMonoidHom)

theorem splitTransportHom_apply (e : Z × U ≃* S) (Φ : A ≃* U × B) (ψ : Z →* B) (s : S) :
    splitTransportHom e Φ ψ s = Φ.symm ((e.symm s).2, ψ (e.symm s).1) := rfl

/-- ★標的側の座標での値 —— これが「`Art` の完全な仕様」である。 -/
theorem apply_splitTransportHom (e : Z × U ≃* S) (Φ : A ≃* U × B) (ψ : Z →* B) (s : S) :
    Φ (splitTransportHom e Φ ψ s) = ((e.symm s).2, ψ (e.symm s).1) := by
  simp [splitTransportHom_apply]

/-- ★**`U` 成分の上での値**（原典の `𝒪^×` 側）。 -/
theorem splitTransportHom_apply_left (e : Z × U ≃* S) (Φ : A ≃* U × B) (ψ : Z →* B) (u : U) :
    splitTransportHom e Φ ψ (e (1, u)) = Φ.symm (u, 1) := by
  simp [splitTransportHom_apply]

/-- ★**`Z` 成分の上での値**（原典の `π` 側）。 -/
theorem splitTransportHom_apply_right (e : Z × U ≃* S) (Φ : A ≃* U × B) (ψ : Z →* B) (z : Z) :
    splitTransportHom e Φ ψ (e (z, 1)) = Φ.symm (1, ψ z) := by
  simp [splitTransportHom_apply]

/-- `ψ` が単射なら `splitTransportHom` も単射。 -/
theorem splitTransportHom_injective (e : Z × U ≃* S) (Φ : A ≃* U × B) {ψ : Z →* B}
    (hψ : Function.Injective ψ) : Function.Injective (splitTransportHom e Φ ψ) := by
  intro s s' h
  rw [splitTransportHom_apply, splitTransportHom_apply] at h
  have h2 := Φ.symm.injective h
  exact e.symm.injective (Prod.ext (hψ (congrArg Prod.snd h2)) (congrArg Prod.fst h2))

/-- ★★**像は `U` 全体 × `ψ` の像**。★第 2 成分が `ψ` の像に潰れているので、
`ψ` が全射でなければ `splitTransportHom` も全射ではない。 -/
theorem range_apply_splitTransportHom (e : Z × U ≃* S) (Φ : A ≃* U × B) (ψ : Z →* B) :
    Set.range (fun s : S => Φ (splitTransportHom e Φ ψ s)) = Set.univ ×ˢ Set.range ψ := by
  ext q
  simp only [Set.mem_range, Set.mem_prod, Set.mem_univ, true_and]
  constructor
  · rintro ⟨s, rfl⟩
    exact ⟨(e.symm s).1, by rw [apply_splitTransportHom]⟩
  · rintro ⟨z, hz⟩
    refine ⟨e (z, q.1), ?_⟩
    rw [apply_splitTransportHom]
    simp [hz]

/-- `Multiplicative ℤ` からの巾写像の像は `{g^n | n : ℤ}`。 -/
theorem range_zpowersHom {G : Type*} [Group G] (g : G) :
    Set.range (zpowersHom G g) = Set.range (fun n : ℤ => g ^ n) := by
  ext x
  constructor
  · rintro ⟨m, rfl⟩; exact ⟨Multiplicative.toAdd m, rfl⟩
  · rintro ⟨n, rfl⟩; exact ⟨Multiplicative.ofAdd n, rfl⟩

/-- ★逆元の巾は同じ集合を掃く（幾何 Frobenius と算術 Frobenius の生成する
部分群が一致すること）。 -/
theorem range_zpow_inv {G : Type*} [Group G] (g : G) :
    Set.range (fun n : ℤ => g⁻¹ ^ n) = Set.range (fun n : ℤ => g ^ n) := by
  ext x
  constructor
  · rintro ⟨n, rfl⟩; exact ⟨-n, by simp [zpow_neg]⟩
  · rintro ⟨n, rfl⟩; exact ⟨-n, by simp [zpow_neg]⟩

/-- ★**抽象核 2**（位相）—— 同相で移した先が稠密なら元も稠密。

★これが `Gal(K^ab/K)` での稠密性に必要な**唯一の**追加である
（`hcont` を落とすと本ファイルでは示せない。§5 参照）。 -/
theorem dense_of_dense_image {A P : Type*} [TopologicalSpace A] [TopologicalSpace P]
    (Φ : A ≃ P) (hcont : Continuous Φ.symm) {s : Set A} (hs : Dense (Φ '' s)) : Dense s := by
  have hcl : Φ.symm '' (closure (Φ '' s)) ⊆ closure (Φ.symm '' (Φ '' s)) :=
    image_closure_subset_closure_image hcont
  rw [hs.closure_eq, Φ.symm_image_image] at hcl
  intro a
  exact hcl ⟨Φ a, Set.mem_univ _, by simp⟩

end AbstractCore

/-! ## §2 一般の Galois 理論 —— 制限写像と 2 つの同型の突き合わせ

★体と Galois の語彙だけ。局所体は出てこない。 -/

section GaloisCore

variable {k E : Type*} [Field k] [Field E] [Algebra k E]

/-- ★**中間体の等式 `A = B` に沿った `Gal` の同一視は、制限写像と可換**。

★#59 対策の要 —— これがあると `Gal(K^ab/K)` の元をすべて
`Γ_K` からの制限（**1 層**）として書ける。 -/
theorem autCongr_equivOfEq_restrictNormalHom
    {A B : IntermediateField k E} (h : A = B) [Normal k A] [Normal k B] (τ : E ≃ₐ[k] E) :
    AlgEquiv.autCongr (IntermediateField.equivOfEq h) (AlgEquiv.restrictNormalHom A τ)
      = AlgEquiv.restrictNormalHom B τ := by
  ext x
  show ((AlgEquiv.restrictNormalHom (A : Type _) τ
      ((IntermediateField.equivOfEq h).symm x) : A) : E)
    = ((AlgEquiv.restrictNormalHom (B : Type _) τ x : B) : E)
  rw [AlgEquiv.restrictNormalHom_apply, AlgEquiv.restrictNormalHom_apply]
  rfl

/-- ★★**`galSupEquivProd` の成分は制限写像である**。

`AbelianDecomposition.lean` の `galSupEquivProd` は第一同型定理を 2 回使って
作られているので、値が「制限の対」であることは別に要る。 -/
theorem galSupEquivProd_restrictNormalHom [IsGalois k E]
    (A B : IntermediateField k E) [Normal k A] [Normal k B] (h : A ⊓ B = ⊥) (τ : E ≃ₐ[k] E) :
    galSupEquivProd A B h (AlgEquiv.restrictNormalHom (A ⊔ B : IntermediateField k E) τ)
      = restrictPairHom A B τ := by
  unfold galSupEquivProd
  simp only [MulEquiv.trans_apply, quotientKerEquivOfSurjective_symm_apply,
    QuotientGroup.quotientMulEquivOfEq_mk]
  rfl

end GaloisCore

/-! ## §3 源の分解 `K^× ≅ π^ℤ × 𝒪_K^×`

原典 §4.2 は `v_L(π_j) = j` と書いて `π_j` を「付値 `j` の標準的な元」として使う。
`K` 上（`L = K`、`π ∈ K` は `ϕ` 不変）ではこれは `π_j = π^j` であり、
`K^× = π^ℤ · 𝒪_K^×` という分解に他ならない。★在庫
`UnitsSplit.lean::unitsSplitEquiv` がそれを与える。 -/

variable {p : ℕ} [Fact p.Prime]

/-- `𝔪 = (π)` なら `π` は既約 —— 離散付値環だから。★これで **Lubin-Tate 級数の
線形係数 `π` と、`K^×` を分裂させる素元が同じもの**になる（原典の `Art_π`）。 -/
theorem irreducible_uniformizer (K : PAdicLocalField p) {ϖ : 𝒪[K.carrier]}
    (h : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {ϖ}) : Irreducible ϖ := by
  haveI := valuationRing_isDVR K
  exact (IsDiscreteValuationRing.irreducible_iff_uniformizer ϖ).mpr h

variable (K : PAdicLocalField p)

/-- **正規化付値** `v : K^× →* ℤ`（`v(ϖ) = 1`）。原典 §4.2 の `v_L`。 -/
noncomputable def artinValuation {ϖ : 𝒪[K.carrier]} (hϖ : Irreducible ϖ) :
    (K.carrier)ˣ →* Multiplicative ℤ :=
  (MonoidHom.fst _ _).comp (unitsSplitEquiv K hϖ).symm.toMonoidHom

/-- **単数部分** `u : K^× →* 𝒪_K^×`（`x = ϖ^{v(x)} · u(x)`）。★`ϖ` に依存する。 -/
noncomputable def artinUnitPart {ϖ : 𝒪[K.carrier]} (hϖ : Irreducible ϖ) :
    (K.carrier)ˣ →* (𝒪[K.carrier])ˣ :=
  (MonoidHom.snd _ _).comp (unitsSplitEquiv K hϖ).symm.toMonoidHom

theorem unitsSplitEquiv_symm_apply {ϖ : 𝒪[K.carrier]} (hϖ : Irreducible ϖ) (x : (K.carrier)ˣ) :
    (unitsSplitEquiv K hϖ).symm x = (artinValuation K hϖ x, artinUnitPart K hϖ x) := rfl

/-- `(1, u) ↦ u`（`𝒪_K^× ⊂ K^×`）。 -/
theorem unitsSplitEquiv_unit {ϖ : 𝒪[K.carrier]} (hϖ : Irreducible ϖ) (u : (𝒪[K.carrier])ˣ) :
    unitsSplitEquiv K hϖ (1, u) = unitsToCarrier K u := by
  apply Units.ext
  simp only [unitsSplitEquiv, MulEquiv.ofBijective_apply, unitsSplitHom, MonoidHom.coe_mk,
    OneHom.coe_mk, Units.coe_map,
    RingHom.toMonoidHom_eq_coe, MonoidHom.coe_coe, toAdd_one, zpow_zero, one_mul]
  rfl

/-- `(1, 1) ↦ ϖ`。 -/
theorem unitsSplitEquiv_gen {ϖ : 𝒪[K.carrier]} (hϖ : Irreducible ϖ) :
    unitsSplitEquiv K hϖ (Multiplicative.ofAdd 1, 1)
      = Units.mk0 (algebraMap 𝒪[K.carrier] K.carrier ϖ)
          (algebraMap_irreducible_ne_zero K hϖ) := by
  apply Units.ext
  simpa [unitsSplitEquiv] using unitsSplitHom_val K hϖ 1 1

/-- ★**退化検査** —— 分解は本物である: `x = ϖ^{v(x)} · u(x)`。 -/
theorem coe_eq_uniformizer_zpow_mul {ϖ : 𝒪[K.carrier]} (hϖ : Irreducible ϖ) (x : (K.carrier)ˣ) :
    (x : K.carrier)
      = (algebraMap 𝒪[K.carrier] K.carrier ϖ) ^ (Multiplicative.toAdd (artinValuation K hϖ x))
        * ((artinUnitPart K hϖ x : 𝒪[K.carrier]) : K.carrier) := by
  have key : unitsSplitEquiv K hϖ (artinValuation K hϖ x, artinUnitPart K hϖ x) = x :=
    MulEquiv.apply_symm_apply _ _
  conv_lhs => rw [← key]
  exact unitsSplitHom_val K hϖ (Multiplicative.toAdd (artinValuation K hϖ x))
    (artinUnitPart K hϖ x)

/-! ## §4 幾何 Frobenius `Frob_K = ϕ^{-1}`

原典 §2.2:「The arithmetic Frobenius ϕ ∈ Gal(K^ur/K) is defined as the element which
reduces mod 𝔭 to the q-th power Frobenius map, and its inverse is denoted by
`Frob_K := ϕ^{-1}` (geometric Frobenius).」 -/

/-- **幾何 Frobenius** `Frob_K = ϕ^{-1}`。★原典 §2.2 の定義そのもの。 -/
noncomputable def geometricFrobenius : unramGal K := (arithFrobenius K)⁻¹

/-- ★**算術 Frobenius は無限位数**（`Gal(K^ur/K) ≅ Ẑ` の非捩れ性のうち、
本ファイルが要る分だけ）。`σ^k` が第 `N` 段を固定するのは `N ∣ k` のときに限る、
という在庫から。 -/
theorem zpow_arithFrobenius_eq_one_iff (k : ℤ) : arithFrobenius K ^ k = 1 ↔ k = 0 := by
  constructor
  · intro h
    by_contra hk
    have hN : k.natAbs + 1 ≠ 0 := Nat.succ_ne_zero _
    have hmem : arithFrobenius K ^ k ∈ (unramLevel K (k.natAbs + 1)).fixingSubgroup := by
      rw [h]; exact one_mem _
    have hdvd := (zpow_arithFrobenius_mem_fixingSubgroup_iff K hN k).mp hmem
    have h2 := Int.natAbs_dvd_natAbs.mpr hdvd
    rw [Int.natAbs_natCast] at h2
    have hpos : 0 < k.natAbs := Int.natAbs_pos.mpr hk
    have := Nat.le_of_dvd hpos h2
    omega
  · rintro rfl; simp

theorem zpow_geometricFrobenius_eq_one_iff (k : ℤ) :
    geometricFrobenius K ^ k = 1 ↔ k = 0 := by
  rw [geometricFrobenius, inv_zpow, inv_eq_one, zpow_arithFrobenius_eq_one_iff]

/-- `ℤ ↪ Gal(K^ur/K)`、`n ↦ Frob_K^n` は単射。 -/
theorem zpowersHom_geometricFrobenius_injective :
    Function.Injective (zpowersHom (unramGal K) (geometricFrobenius K)) := by
  intro m n h
  have h1 : geometricFrobenius K ^ (Multiplicative.toAdd m)
      = geometricFrobenius K ^ (Multiplicative.toAdd n) := h
  have h2 : geometricFrobenius K ^ (Multiplicative.toAdd m - Multiplicative.toAdd n) = 1 := by
    rw [zpow_sub, h1, mul_inv_cancel]
  have h3 := (zpow_geometricFrobenius_eq_one_iff K _).mp h2
  have h4 : Multiplicative.toAdd m = Multiplicative.toAdd n := by omega
  exact Multiplicative.toAdd.injective h4

/-! ## §5 標的の分解 `Gal(M/K) ≅ 𝒪_K^× × Gal(K^ur/K)`

★`M` は「`K_π · K^ur` に等しい中間体」としてパラメータ化してある。
§6 で `M := K^ab` を **LKW（Theorem 6.15）** で代入する。 -/

variable [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {M : IntermediateField K.carrier K.closure}

/-- ★★**`Gal(M/K) ≅ 𝒪_K^× × Gal(K^ur/K)`**（`M = K_π · K^ur` のとき）。

`AbelianDecomposition.lean::lubinTateUnramifiedGalEquivProd` を
中間体の等式 `hM` に沿って移しただけ。 -/
noncomputable def abelianGalEquivProd
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) :
    (M ≃ₐ[K.carrier] M) ≃* (𝒪[K.carrier])ˣ × unramGal K :=
  (AlgEquiv.autCongr (IntermediateField.equivOfEq hM)).symm.trans
    (lubinTateUnramifiedGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf)

/-- ★★★**`Φ` の成分は `Γ_K` からの制限そのもの**。

    Φ (τ|_M) = ( ρ(τ|_{K_π}), τ|_{K^ur} )

★#59 回避（定型 (b)）: 左辺・右辺のすべてが `K̄` からの **1 層の**制限で書かれており、
中間体の中の中間体は現れない。 -/
theorem abelianGalEquivProd_restrictNormalHom [Normal K.carrier M]
    [Normal K.carrier (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf)]
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) (τ : K.absGal) :
    abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM
        (AlgEquiv.restrictNormalHom (M : Type _) τ)
      = (lubinTateClosureGalEquivUnits K hq hπmax hπne0 f hf0 hf1 hf
            (AlgEquiv.restrictNormalHom
              (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf : Type _) τ),
          AlgEquiv.restrictNormalHom (unramifiedClosure K : Type _) τ) := by
  haveI := isGalois_carrier_closure K
  haveI : Normal K.carrier
      (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
        IntermediateField K.carrier K.closure) := hM ▸ (inferInstance : Normal K.carrier M)
  have h1 : (AlgEquiv.autCongr (IntermediateField.equivOfEq hM)).symm
      (AlgEquiv.restrictNormalHom (M : Type _) τ)
      = AlgEquiv.restrictNormalHom
          (((lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
            IntermediateField K.carrier K.closure) : Type _)) τ := by
    rw [← autCongr_equivOfEq_restrictNormalHom hM τ]
    exact MulEquiv.symm_apply_apply _ _
  unfold abelianGalEquivProd lubinTateUnramifiedGalEquivProd
  rw [MulEquiv.trans_apply, h1, MulEquiv.trans_apply, galSupEquivProd_restrictNormalHom]
  rfl

/-! ### `Φ` は位相群の同型でもある

★これがあると像の稠密性が `Gal(M/K)` そのものの中で言える（§7）。 -/

/-- ★★**`Φ` は連続**。

`Γ_K ↠ Gal(M/K)` はコンパクト → Hausdorff の連続全射なので**商写像**であり、
`Φ` の連続性は `Φ ∘ (制限)` の連続性に帰着する。後者は
`abelianGalEquivProd_restrictNormalHom` により
`(ρ ∘ 制限, 制限)` に等しく、両成分とも連続
（`InfiniteGalois.restrictNormalHom_continuous` と
`LubinTateClosureTopology.lean::continuous_lubinTateClosureGalEquivUnits`）。 -/
theorem continuous_abelianGalEquivProd [IsGalois K.carrier M]
    [Normal K.carrier (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf)]
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) :
    Continuous (abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM) := by
  haveI := isGalois_carrier_closure K
  haveI := normal_unramifiedClosure K
  have hres : Continuous (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      (M : Type _)) := InfiniteGalois.restrictNormalHom_continuous M
  have hsurj := AlgEquiv.restrictNormalHom_surjective (F := K.carrier)
    (K₁ := (M : Type _)) K.closure
  have hquot : Topology.IsQuotientMap (AlgEquiv.restrictNormalHom (F := K.carrier)
      (K₁ := K.closure) (M : Type _)) := (hres.isClosedMap).isQuotientMap hres hsurj
  rw [hquot.continuous_iff]
  have heq : (fun τ : K.absGal => abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM
        (AlgEquiv.restrictNormalHom (M : Type _) τ))
      = fun τ : K.absGal =>
        (lubinTateClosureGalEquivUnits K hq hπmax hπne0 f hf0 hf1 hf
            (AlgEquiv.restrictNormalHom
              (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf : Type _) τ),
          AlgEquiv.restrictNormalHom (unramifiedClosure K : Type _) τ) :=
    funext (abelianGalEquivProd_restrictNormalHom K hq hπmax hπne0 f hf0 hf1 hf hM)
  show Continuous (fun τ : K.absGal => abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM
    (AlgEquiv.restrictNormalHom (M : Type _) τ))
  rw [heq]
  exact ((continuous_lubinTateClosureGalEquivUnits K hq hπmax hπne0 f hf0 hf1 hf).comp
      (InfiniteGalois.restrictNormalHom_continuous _)).prodMk
    (InfiniteGalois.restrictNormalHom_continuous _)

/-- ★★**`Φ.symm` も連続** —— `Gal(M/K)` はコンパクト、`𝒪_K^× × Gal(K^ur/K)` は
Hausdorff なので `Continuous.homeoOfEquivCompactToT2` が逆写像の連続性をくれる。
★これで `Φ` は**位相群の同型**である。 -/
theorem continuous_abelianGalEquivProd_symm [IsGalois K.carrier M]
    [Normal K.carrier (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf)]
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) :
    Continuous (abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM).symm := by
  have hcont : Continuous
      ⇑(abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM).toEquiv :=
    continuous_abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM
  exact hcont.homeoOfEquivCompactToT2.symm.continuous

/-! ## §6 Artin 写像 `Art_π` -/

def artinMap.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Definition 4.10", sectionId := "def-4-10" }

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**Artin 写像 `Art_π : K^× →* Gal(M/K)`**（`M = K_π · K^ur`）。

原文 (Yoshida08 p.9):
> Definition 4.10. For any f ∈ O[scr]_L[X] with L/K finite, set K^m := K^urL^m_f. Then K^m/K is finitely ramified, and Galois by Proposition 4.7(i). By Lemma 2.2, the completion of K^m is KL^m_f = K^m and K^m = K^m ∩ K^sep, thus independent of f. Setting K^LT := _m≥1 K^m = K^LT ∩K^sep, we have W(K^LT/K) ∼ = W(K^LT/K) by the remark after Definition 2.5. We call a finite extension of K a Lubin-Tate extension if it is contained in K^LT. We call the inverse of ρ the Artin map of K and write Art_K : K^× ∼ =−→ W(K^LT/K). We have v ◦ Art_K = v.

定義は 2 つの直積分解の合成:

    K^×  ≅  π^ℤ × 𝒪_K^×  ∋ (π^n, u)  ↦  (u, Frob_K^n) ∈ 𝒪_K^× × Gal(K^ur/K)  ≅  Gal(M/K)

★★**`π` に依存する**（引数 `hπmax`・`f`・`hf0`・`hf1`・`hf` に現れている）。
★**非依存性は主張していない**（Corollary 4.9 + Dwork を要する別の節点）。 -/
noncomputable def artinMap
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) :
    (K.carrier)ˣ →* (M ≃ₐ[K.carrier] M) :=
  splitTransportHom (unitsSplitEquiv K (irreducible_uniformizer K hπmax))
    (abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM)
    (zpowersHom (unramGal K) (geometricFrobenius K))

/-- ★★**`Art_π` の完全な仕様**（標的の座標での値）:

    Φ (Art_π x) = ( u(x), Frob_K^{v(x)} ). -/
theorem abelianGalEquivProd_artinMap
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) (x : (K.carrier)ˣ) :
    abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM
        (artinMap K hq hπmax hπne0 f hf0 hf1 hf hM x)
      = (artinUnitPart K (irreducible_uniformizer K hπmax) x,
          geometricFrobenius K ^
            (Multiplicative.toAdd (artinValuation K (irreducible_uniformizer K hπmax) x))) :=
  apply_splitTransportHom _ _ _ x

/-- ★★**原典 Definition 4.10 の `v ∘ Art_K = v`** ——
`Art_π(x)` の `K^ur` 成分はちょうど `Frob_K^{v(x)}`（幾何 Frobenius の `v(x)` 乗）。 -/
theorem artinMap_unramified_component
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) (x : (K.carrier)ˣ) :
    (abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM
        (artinMap K hq hπmax hπne0 f hf0 hf1 hf hM x)).2
      = geometricFrobenius K ^
          (Multiplicative.toAdd (artinValuation K (irreducible_uniformizer K hπmax) x)) := by
  rw [abelianGalEquivProd_artinMap]

/-- ★★**`𝒪_K^×` の上での値** —— `Art_π(u)` は `K_π` 上で相互律同型の逆、
`K^ur` 上で恒等。★原典 Proposition 4.7(ii) の `j = 0` の場合。 -/
theorem artinMap_unitsToCarrier
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) (u : (𝒪[K.carrier])ˣ) :
    artinMap K hq hπmax hπne0 f hf0 hf1 hf hM (unitsToCarrier K u)
      = (abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM).symm (u, 1) := by
  rw [← unitsSplitEquiv_unit K (irreducible_uniformizer K hπmax) u]
  exact splitTransportHom_apply_left _ _ _ u

/-- ★★**素元 `π` の上での値** —— `Art_π(π)` は `K_π` 上で恒等、
`K^ur` 上で幾何 Frobenius `Frob_K = ϕ^{-1}`。
★原典 Proposition 4.7(ii) の `v(x) = 1`（したがって `j = −1`、`σ|_{K̂} = ϕ^{-1}`）。 -/
theorem artinMap_uniformizer
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) :
    artinMap K hq hπmax hπne0 f hf0 hf1 hf hM
        (Units.mk0 (algebraMap 𝒪[K.carrier] K.carrier π)
          (algebraMap_irreducible_ne_zero K (irreducible_uniformizer K hπmax)))
      = (abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM).symm
          (1, geometricFrobenius K) := by
  rw [← unitsSplitEquiv_gen K (irreducible_uniformizer K hπmax)]
  have h := splitTransportHom_apply_right (unitsSplitEquiv K (irreducible_uniformizer K hπmax))
    (abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM)
    (zpowersHom (unramGal K) (geometricFrobenius K)) (Multiplicative.ofAdd 1)
  rw [show artinMap K hq hπmax hπne0 f hf0 hf1 hf hM = splitTransportHom
      (unitsSplitEquiv K (irreducible_uniformizer K hπmax))
      (abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM)
      (zpowersHom (unramGal K) (geometricFrobenius K)) from rfl, h]
  congr 1

/-- ★★**`Art_π` は単射**（原典の `Art_K : K^× ≅ W(K^LT/K)` の単射性）。 -/
theorem artinMap_injective
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) :
    Function.Injective (artinMap K hq hπmax hπne0 f hf0 hf1 hf hM) :=
  splitTransportHom_injective _ _ (zpowersHom_geometricFrobenius_injective K)

/-- ★★★**`Art_π` の特徴づけ**（`Γ_K` の言葉で）——
`τ|_M = Art_π(x)` であるための必要十分条件は、
`τ` が `K_π` 上で `u(x)` 倍作用し、`K^ur` 上で `Frob_K^{v(x)}` であること。

★**これが「`𝒪^×` 上と `π` 上の両方を決めた」ことの内容**である。 -/
theorem restrictNormalHom_eq_artinMap_iff [Normal K.carrier M]
    [Normal K.carrier (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf)]
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) (x : (K.carrier)ˣ) (τ : K.absGal) :
    AlgEquiv.restrictNormalHom (M : Type _) τ = artinMap K hq hπmax hπne0 f hf0 hf1 hf hM x
      ↔ (lubinTateClosureGalEquivUnits K hq hπmax hπne0 f hf0 hf1 hf
              (AlgEquiv.restrictNormalHom
                (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf : Type _) τ)
            = artinUnitPart K (irreducible_uniformizer K hπmax) x
          ∧ AlgEquiv.restrictNormalHom (unramifiedClosure K : Type _) τ
            = geometricFrobenius K ^
                (Multiplicative.toAdd (artinValuation K (irreducible_uniformizer K hπmax) x))) := by
  rw [← (abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM).injective.eq_iff,
    abelianGalEquivProd_restrictNormalHom, abelianGalEquivProd_artinMap, Prod.ext_iff]

/-! ## §7 像 —— Weil 群であって `Gal` 全体ではない -/

/-- ★★**像は `𝒪_K^× × Frob^ℤ`**。★第 2 成分が `Frob^ℤ` に留まる
（`Gal(K^ur/K) ≅ Ẑ` の全体ではない）ことが、`Art_π` が**全射でない**理由である。 -/
theorem range_abelianGalEquivProd_artinMap
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) :
    Set.range (fun x : (K.carrier)ˣ => abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM
        (artinMap K hq hπmax hπne0 f hf0 hf1 hf hM x))
      = Set.univ ×ˢ Set.range (fun n : ℤ => arithFrobenius K ^ n) := by
  rw [artinMap, range_apply_splitTransportHom, range_zpowersHom, geometricFrobenius,
    range_zpow_inv]

/-- ★★★**像はちょうど Weil 群** `W(M/K) = {σ | σ|_{K^ur} ∈ Frob^ℤ}`。

原典 Definition 4.10 の `Art_K : K^× ≅ W(K^LT/K)` の「像の側」。
★**`Gal(M/K)` 全体ではない**。 -/
theorem range_artinMap
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) :
    Set.range (artinMap K hq hπmax hπne0 f hf0 hf1 hf hM)
      = {σ : M ≃ₐ[K.carrier] M | ∃ n : ℤ,
          (abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM σ).2
            = arithFrobenius K ^ n} := by
  ext σ
  constructor
  · rintro ⟨x, rfl⟩
    have hmem : abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM
        (artinMap K hq hπmax hπne0 f hf0 hf1 hf hM x)
        ∈ Set.range (fun x : (K.carrier)ˣ =>
          abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM
            (artinMap K hq hπmax hπne0 f hf0 hf1 hf hM x)) := ⟨x, rfl⟩
    rw [range_abelianGalEquivProd_artinMap] at hmem
    obtain ⟨n, hn⟩ := hmem.2
    exact ⟨n, hn.symm⟩
  · rintro ⟨n, hn⟩
    have hmem : abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM σ
        ∈ Set.univ ×ˢ Set.range (fun n : ℤ => arithFrobenius K ^ n) :=
      ⟨Set.mem_univ _, ⟨n, hn.symm⟩⟩
    rw [← range_abelianGalEquivProd_artinMap K hq hπmax hπne0 f hf0 hf1 hf hM] at hmem
    obtain ⟨x, hx⟩ := hmem
    exact ⟨x, (abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM).injective hx⟩

/-- ★★**像は `Φ` で移した先で稠密**。

`𝒪_K^×` 成分は全体、`Gal(K^ur/K)` 成分は `Frob^ℤ` で、後者は
`UnramifiedZhat.lean::dense_zpowers_frobenius` により稠密。 -/
theorem dense_range_abelianGalEquivProd_artinMap
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) :
    Dense (Set.range (fun x : (K.carrier)ˣ =>
      abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM
        (artinMap K hq hπmax hπne0 f hf0 hf1 hf hM x))) := by
  rw [range_abelianGalEquivProd_artinMap]
  exact dense_univ.prod (dense_zpowers_frobenius K (arithFrobenius_mem_unramLevelGeneratorSet K))

/-- ★★★**`Gal(M/K)` そのものの中で像は稠密**（仮定なし）。

★`Φ` が位相群の同型であること（`continuous_abelianGalEquivProd_symm`）を使って
`Φ` で移した先の稠密性を引き戻すだけ。
★★**稠密であって全射ではない** —— `range_artinMap` が像を Weil 群と同定している。 -/
theorem dense_range_artinMap [IsGalois K.carrier M]
    [Normal K.carrier (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf)]
    (hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = M) :
    Dense (Set.range (artinMap K hq hπmax hπne0 f hf0 hf1 hf hM)) := by
  refine dense_of_dense_image
    (abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM).toEquiv
    (continuous_abelianGalEquivProd_symm K hq hπmax hπne0 f hf0 hf1 hf hM) ?_
  rw [← Set.range_comp]
  exact dense_range_abelianGalEquivProd_artinMap K hq hπmax hπne0 f hf0 hf1 hf hM

/-! ## §8 `K^ab` への特殊化 —— ★ここで LKW（Theorem 6.15）を消費する -/

/-- ★★★**`Art_π : K^× →* Gal(K^ab/K)`**。

★**LKW を消費している**: `abelianClosure_eq_lubinTateClosure_sup_unramifiedClosure`
（Theorem 6.15、`LocalClassFieldTheory.lean`）が `K^ab = K_π · K^ur` を与え、
それを `artinMap` の `hM` に代入している。★これが無いと標的は `Gal(K_π·K^ur/K)` 止まり。 -/
noncomputable def artinMapAbelian :
    (K.carrier)ˣ →* (abelianClosure K ≃ₐ[K.carrier] abelianClosure K) :=
  artinMap K hq hπmax hπne0 f hf0 hf1 hf
    (abelianClosure_eq_lubinTateClosure_sup_unramifiedClosure K hq hπmax hπne0 f hf0 hf1 hf).symm

theorem artinMapAbelian_injective : Function.Injective
    (artinMapAbelian K hq hπmax hπne0 f hf0 hf1 hf) :=
  artinMap_injective K hq hπmax hπne0 f hf0 hf1 hf
    (abelianClosure_eq_lubinTateClosure_sup_unramifiedClosure K hq hπmax hπne0 f hf0 hf1 hf).symm

def exists_artinMap.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Definition 4.10", sectionId := "def-4-10" }

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**Yoshida 2008 Definition 4.10 —— 仮説なしの形の Artin 写像。**

原文 (Yoshida08 p.9):
> Definition 4.10. For any f ∈ O[scr]_L[X] with L/K finite, set K^m := K^urL^m_f. Then K^m/K is finitely ramified, and Galois by Proposition 4.7(i). By Lemma 2.2, the completion of K^m is KL^m_f = K^m and K^m = K^m ∩ K^sep, thus independent of f. Setting K^LT := _m≥1 K^m = K^LT ∩K^sep, we have W(K^LT/K) ∼ = W(K^LT/K) by the remark after Definition 2.5. We call a finite extension of K a Lubin-Tate extension if it is contained in K^LT. We call the inverse of ρ the Artin map of K and write Art_K : K^× ∼ =−→ W(K^LT/K). We have v ◦ Art_K = v.

任意の p 進局所体 `K` について、同型 `Φ : Gal(K^ab/K) ≅ 𝒪_K^× × Gal(K^ur/K)` と
単射準同型 `Art : K^× →* Gal(K^ab/K)` と素元 `ϖ ∈ K^×` が在って

* `Art(u) = Φ^{-1}(u, 1)`（`u ∈ 𝒪_K^×`）
* `Art(ϖ) = Φ^{-1}(1, Frob_K)`（`Frob_K` は**幾何** Frobenius）
* 像は `Φ` で移すと稠密

★`π`・`f` の選択は `∃` の内側に閉じ込めてある（`Art` は選択に依存する。
★**依存しないとは主張していない**）。 -/
theorem exists_artinMap (K : PAdicLocalField p) :
    ∃ (Φ : (abelianClosure K ≃ₐ[K.carrier] abelianClosure K) ≃*
        ((𝒪[K.carrier])ˣ × unramGal K))
      (Art : (K.carrier)ˣ →* (abelianClosure K ≃ₐ[K.carrier] abelianClosure K))
      (ϖ : (K.carrier)ˣ),
      Function.Injective Art
      ∧ (∀ u : (𝒪[K.carrier])ˣ, Art (unitsToCarrier K u) = Φ.symm (u, 1))
      ∧ Art ϖ = Φ.symm (1, geometricFrobenius K)
      ∧ Dense (Set.range fun x : (K.carrier)ˣ => Φ (Art x))
      ∧ Dense (Set.range Art) := by
  haveI := isAdicComplete_valuationRing K
  haveI := valuationRing_isDVR K
  obtain ⟨ϖ, hϖirr⟩ := IsDiscreteValuationRing.exists_irreducible (𝒪[K.carrier])
  have hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {ϖ} :=
    (IsDiscreteValuationRing.irreducible_iff_uniformizer ϖ).mp hϖirr
  have hπne0 : ϖ ≠ 0 := hϖirr.ne_zero
  have hq : Fintype.card 𝓀[K.carrier] = p ^ (absoluteInertiaDegree K) := by
    rw [← Nat.card_eq_fintype_card]
    exact residueCard_eq_pow K
  obtain ⟨f, hf0, hf1, hf⟩ := exists_lubinTateSeries (A := 𝒪[K.carrier]) hq hπmax
  haveI := normal_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf
  have hM : (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) = abelianClosure K :=
    (abelianClosure_eq_lubinTateClosure_sup_unramifiedClosure K hq hπmax hπne0 f hf0 hf1 hf).symm
  refine ⟨abelianGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf hM,
    artinMap K hq hπmax hπne0 f hf0 hf1 hf hM,
    Units.mk0 (algebraMap 𝒪[K.carrier] K.carrier ϖ)
      (algebraMap_irreducible_ne_zero K (irreducible_uniformizer K hπmax)),
    artinMap_injective K hq hπmax hπne0 f hf0 hf1 hf hM,
    fun u => artinMap_unitsToCarrier K hq hπmax hπne0 f hf0 hf1 hf hM u,
    artinMap_uniformizer K hq hπmax hπne0 f hf0 hf1 hf hM,
    dense_range_abelianGalEquivProd_artinMap K hq hπmax hπne0 f hf0 hf1 hf hM,
    dense_range_artinMap K hq hπmax hπne0 f hf0 hf1 hf hM⟩

end ABC3.Found.PGC
