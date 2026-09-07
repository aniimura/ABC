import ABC3.Found.PGC.Prop22FilteredHypothesis
import ABC3.Found.PGC.RamificationImageStage

/-!
# 分岐フィルトレーションの**自然性** —— 段データ 1 段に還元する

`Found/PGC/RamificationNaturality.lean::IsNaturalFiltration` は

```
∀ {K K'} (α : K.carrier ≃ₐ[ℚ_p] K'.carrier) (v : ℝ),
  Subgroup.map (galContinuousMulEquiv α).toMulEquiv (RF.Gv K v) = RF.Gv K' v
```

である。★自由な `RF` では**偽**(D32、退化した `Gv ≡ ⊤` が絡む)なので、
本ファイルは実物 `ABC3.Found.PGC.ramificationFiltration p`
(`Found/PGC/UnramifiedBaseChangeInvariance.lean`、仮定ゼロ)に固定して扱う。

## ★★何が入ったか

### §1 抽象核 —— `PAdicLocalField`・分岐・付値・Galois が **1 語も出ない**

* `map_eq_map_of_conj`(**純群論**) —— `H` が正規なら、
  「共役でずれるだけの 2 つの同型」`Φ, Ψ` の像は**同じ**。
* `map_iInf_eq_of_bij`(**純群論**) —— 同型 `Φ` が添字集合 `S ↔ S'` を
  対応させ、各段で `map Φ (A N) = A' (map Φ N)` なら、**逆極限も対応する**。
* `map_mem_openNormalBase` / `comap_mem_openNormalBase`(**位相群のみ**) ——
  位相群の同型は開正規部分群を開正規部分群に写す。

### §2 ★★延長 `ᾱ` の選択に依らない(どの `RamificationFiltration` でも)

`map_gv_eq_of_ext` —— `α` の**任意の**互換な延長 `ᾱ` から作った同型は
`galContinuousMulEquiv α` と**同じ像**を与える。中身は §1 の `map_eq_map_of_conj`
＋ `galMulEquivOf_indep`(内部自己同型のずれしか生まない)＋ `RF.isNormal` だけ。

★これで「`IsNaturalFiltration` を示すには**都合のよい延長を 1 つ選べばよい**」
ことが分かる。★★**原典が `Out`(外部同型)で述べている理由がここに出ている。**

### §3 ★★★義務を段データ 1 段に落とす

`IsNaturalStage p` —— 「開正規部分群 `N` と `v > 0` に対し
`map Φ (absGalStage K N v) = absGalStage K' (map Φ N) v`」。

`isNaturalFiltration_of_isNaturalStage : IsNaturalStage p →
IsNaturalFiltration (ramificationFiltration p)`。

★逆も真(`isNaturalStage_of_isNaturalFiltration` は**取っていない**——
`Γ^v · M = S M v` から段を復元する向きは `compat` を要するので、本ファイルでは
片側だけを主張する)。

★直前の波の `isNaturalFiltration_ramificationFiltration_iff_pos`
(`Found/PGC/Prop22FilteredHypothesis.lean`)で `v ≤ 0` は既に落ちているので、
**`v > 0` だけの仮説**にしてある。

### §4 ★★非空虚性(無条件) —— 仮説は空ではない

* `map_absGalStage_of_nonpos` —— ★**`v ≤ 0` では段データの自然性は無条件に成り立つ**
  (`absGalStage K N v = I_K ⊔ N` と Corollary 1.3)。★`v > 0` に絞ったことが正しいことの裏。
* `map_absGalStage_inner` / `map_gv_inner` —— ★**内部自己同型では無条件**
  (`absGalStage K N v` が正規で `N` も正規だから)。
* `map_gv_refl` —— ★★**`α = id` の場合は無条件に成り立つ**。
  ★これは自明ではない: `galContinuousMulEquiv (AlgEquiv.refl)` は
  `extendToClosure` の選択のせいで**恒等写像ではない**(`Γ_K` の内部自己同型である)。

## ★★★残った穴 —— ★**同じシフトの後続ファイルで埋まった**

★★**`IsNaturalStage p` は `Found/PGC/StageTransport.lean` で有限段の主張
`IsNaturalStageUpper p` に落とし、`Found/PGC/StageUpperNaturality.lean` で
**無条件に証明した**(`isNaturalFiltration_ramificationFiltration`)。
★本ファイルの `IsNaturalStage` は「途中の仮説」として残してある(依存の向きの都合)。

以下は本ファイルを書いた時点での見通しであり、実際にその通りになった:

残っているのは
「体の同型 `α` の延長 `ᾱ` が有限段 `L = K(x)` の**上付き分岐群**を
`L' = K'(ᾱ x)` の上付き分岐群に写す」という**具体層**であり、
`Found/PGC/UnramifiedBaseChangeInvariance.lean::map_upperRamificationGroup_eq`
(`ramIndex` が対応すれば `G^m` も対応する)に還元されるが、そのためには

1. `ᾱ` が `𝒪_L ≃+* 𝒪_{L'}` を誘導すること(`norm_extendToClosure` から出るはず)、
2. その環同型が `Gal(L/K) ≃ Gal(L'/K')` と同変であること、
3. `addVal` が保たれること

の 3 本が要る。★**本ファイルでは着手していない(測っていない、ではなく着手していない)。**

## 逸脱の記録

1. `Skeleton/**` は 1 行も書き換えていない。
2. `IsNaturalStage` は `v > 0` にだけ条件を課す。`v ≤ 0` が自動であることは
   `map_absGalStage_of_nonpos` で**証明した**ので、これは弱め方ではなく等価な形である
   (`isNaturalStage_iff_forall` を参照)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC

/-! ## §1 ★★抽象核 —— 分岐・付値・Galois・`ℚ_p` が 1 語も出てこない -/

section AbstractCore

variable {Γ Γ' : Type*} [Group Γ] [Group Γ']

/-- ★★**抽象核(純群論)** —— `H` が**正規**なら、「共役でずれるだけの 2 つの同型」
`Φ`, `Ψ` は `H` の像として**同じもの**を与える。

★これが「外部同型 `Out` で述べれば十分」の中身である。 -/
theorem map_eq_map_of_conj (Φ Ψ : Γ ≃* Γ') (H : Subgroup Γ) [H.Normal]
    (h : ∀ g : Γ, ∃ c : Γ', Φ g = c * Ψ g * c⁻¹) :
    Subgroup.map (Φ : Γ →* Γ') H = Subgroup.map (Ψ : Γ →* Γ') H := by
  haveI hΨ : (Subgroup.map (Ψ : Γ →* Γ') H).Normal :=
    Subgroup.Normal.map inferInstance _ Ψ.surjective
  haveI hΦ : (Subgroup.map (Φ : Γ →* Γ') H).Normal :=
    Subgroup.Normal.map inferInstance _ Φ.surjective
  apply le_antisymm
  · rintro _ ⟨x, hx, rfl⟩
    obtain ⟨c, hc⟩ := h x
    show Φ x ∈ Subgroup.map (Ψ : Γ →* Γ') H
    rw [hc]
    exact hΨ.conj_mem _ ⟨x, hx, rfl⟩ c
  · rintro _ ⟨x, hx, rfl⟩
    obtain ⟨c, hc⟩ := h x
    show Ψ x ∈ Subgroup.map (Φ : Γ →* Γ') H
    have hc' : Ψ x = c⁻¹ * Φ x * c := by rw [hc]; group
    rw [hc']
    simpa using hΦ.conj_mem _ ⟨x, hx, rfl⟩ c⁻¹

/-- ★★**抽象核(純群論)** —— 同型による像は**逆極限と可換**。

添字集合 `S`・`S'` が `Φ` で対応し(`hmap` / `hcomap`)、各段で
`map Φ (A N) = A' (map Φ N)` なら、`⨅` も対応する。 -/
theorem map_iInf_eq_of_bij (Φ : Γ ≃* Γ') {S : Set (Subgroup Γ)} {S' : Set (Subgroup Γ')}
    (A : Subgroup Γ → Subgroup Γ) (A' : Subgroup Γ' → Subgroup Γ')
    (hmap : ∀ ⦃N⦄, N ∈ S → Subgroup.map (Φ : Γ →* Γ') N ∈ S')
    (hcomap : ∀ ⦃N'⦄, N' ∈ S' → Subgroup.comap (Φ : Γ →* Γ') N' ∈ S)
    (hA : ∀ ⦃N⦄, N ∈ S → Subgroup.map (Φ : Γ →* Γ') (A N) = A' (Subgroup.map (Φ : Γ →* Γ') N)) :
    Subgroup.map (Φ : Γ →* Γ') (⨅ N : S, A (N : Subgroup Γ))
      = ⨅ N' : S', A' (N' : Subgroup Γ') := by
  have hcm : ∀ N' : Subgroup Γ',
      Subgroup.map (Φ : Γ →* Γ') (Subgroup.comap (Φ : Γ →* Γ') N') = N' :=
    fun N' => Subgroup.map_comap_eq_self_of_surjective Φ.surjective N'
  apply le_antisymm
  · rintro _ ⟨x, hx, rfl⟩
    simp only [SetLike.mem_coe, Subgroup.mem_iInf, Subtype.forall] at hx
    simp only [Subgroup.mem_iInf, Subtype.forall]
    intro N' hN'
    rw [← hcm N', ← hA (hcomap hN')]
    exact ⟨x, hx _ (hcomap hN'), rfl⟩
  · intro y hy
    simp only [Subgroup.mem_iInf, Subtype.forall] at hy
    obtain ⟨x, rfl⟩ := Φ.surjective y
    refine ⟨x, ?_, rfl⟩
    simp only [SetLike.mem_coe, Subgroup.mem_iInf, Subtype.forall]
    intro N hN
    have h1 := hy _ (hmap hN)
    rw [← hA hN] at h1
    obtain ⟨z, hz, hzeq⟩ := h1
    rwa [Φ.injective hzeq] at hz

end AbstractCore

section AbstractTopo

variable {Γ Γ' : Type*} [Group Γ] [TopologicalSpace Γ] [Group Γ'] [TopologicalSpace Γ']

/-- ★**抽象核(位相群のみ)** —— 位相群の同型は開正規部分群を開正規部分群に写す。 -/
theorem map_mem_openNormalBase (Φ : ContinuousMulEquiv Γ Γ') {N : Subgroup Γ}
    (hN : N ∈ openNormalBase Γ) : Subgroup.map (Φ.toMulEquiv : Γ →* Γ') N ∈ openNormalBase Γ' := by
  refine ⟨?_, Subgroup.Normal.map hN.2 _ Φ.surjective⟩
  rw [Subgroup.coe_map]
  exact Φ.toHomeomorph.isOpenMap _ hN.1

/-- ★**抽象核(位相群のみ)** —— 引き戻しも同じ。 -/
theorem comap_mem_openNormalBase (Φ : ContinuousMulEquiv Γ Γ') {N' : Subgroup Γ'}
    (hN' : N' ∈ openNormalBase Γ') :
    Subgroup.comap (Φ.toMulEquiv : Γ →* Γ') N' ∈ openNormalBase Γ :=
  ⟨hN'.1.preimage Φ.continuous_toFun, hN'.2.comap _⟩

end AbstractTopo

variable {p : ℕ} [Fact p.Prime]

/-! ## §2 ★★延長 `ᾱ` の選択に依らない —— どの `RamificationFiltration` でも -/

/-- ★★**`Γ_K^v` の像は延長 `ᾱ` の取り方に依らない**。

`galMulEquivOf_indep`(2 つの延長は内部自己同型でずれるだけ)と
`RF.isNormal`(`Γ_K^v` は正規)を §1 の `map_eq_map_of_conj` に通しただけ。

★★これで `IsNaturalFiltration` を示すときに**都合のよい延長を 1 つ選んでよい**。 -/
theorem map_gv_eq_of_ext (RF : RamificationFiltration p) {K K' : PAdicLocalField p}
    (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (ᾱ : K.closure ≃+* K'.closure)
    (hfwd : ∀ x : K.carrier, ᾱ (algebraMap K.carrier K.closure x)
      = algebraMap K'.carrier K'.closure (α x)) (v : ℝ) :
    Subgroup.map ((galMulEquivOf α ᾱ hfwd : K.absGal ≃* K'.absGal) : K.absGal →* K'.absGal)
        (RF.Gv K v)
      = Subgroup.map (((galContinuousMulEquiv α).toMulEquiv : K.absGal ≃* K'.absGal) :
          K.absGal →* K'.absGal) (RF.Gv K v) := by
  haveI := RF.isNormal K v
  exact map_eq_map_of_conj _ _ _ (galMulEquiv_conj_indep α ᾱ hfwd)

def map_gv_eq_of_ext.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

/-! ## §3 ★★★義務を段データ 1 段に落とす -/

/-- 実物の分岐フィルトレーションは段データの**逆極限**である(1 層の展開、`rfl`)。 -/
theorem ramificationFiltration_Gv_eq_iInf (K : PAdicLocalField p) (v : ℝ) :
    (ramificationFiltration p).Gv K v
      = ⨅ N : (openNormalBase K.absGal), absGalStage K (N : Subgroup K.absGal) v := rfl

/-- ★★★**段データの自然性** —— `IsNaturalFiltration (ramificationFiltration p)` の
**1 段小さい**仮説。

`v > 0` にしか条件を課さない(`v ≤ 0` は `map_absGalStage_of_nonpos` で無条件に
成り立つ)。 -/
def IsNaturalStage (p : ℕ) [Fact p.Prime] : Prop :=
  ∀ {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier)
    ⦃N : Subgroup K.absGal⦄, N ∈ openNormalBase K.absGal → ∀ v : ℝ, 0 < v →
      Subgroup.map (((galContinuousMulEquiv α).toMulEquiv : K.absGal ≃* K'.absGal) :
            K.absGal →* K'.absGal) (absGalStage K N v)
        = absGalStage K' (Subgroup.map (((galContinuousMulEquiv α).toMulEquiv :
            K.absGal ≃* K'.absGal) : K.absGal →* K'.absGal) N) v

def IsNaturalStage.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

/-- ★★★★**段データの自然性から分岐フィルトレーションの自然性が出る**。

中身は §1 の `map_iInf_eq_of_bij` ＋ `map_mem_openNormalBase` ＋
直前の波の `isNaturalFiltration_ramificationFiltration_iff_pos`(`v ≤ 0` を落とす)。 -/
theorem isNaturalFiltration_of_isNaturalStage (h : IsNaturalStage p) :
    IsNaturalFiltration (ramificationFiltration p) := by
  refine isNaturalFiltration_ramificationFiltration_iff_pos.mpr ?_
  intro K K' α v hv
  have key := map_iInf_eq_of_bij (Γ := K.absGal) (Γ' := K'.absGal)
    ((galContinuousMulEquiv α).toMulEquiv)
    (S := openNormalBase K.absGal) (S' := openNormalBase K'.absGal)
    (fun N => absGalStage K N v) (fun N' => absGalStage K' N' v)
    (fun _ hN => map_mem_openNormalBase (galContinuousMulEquiv α) hN)
    (fun _ hN' => comap_mem_openNormalBase (galContinuousMulEquiv α) hN')
    (fun _ hN => h α hN v hv)
  rw [ramificationFiltration_Gv_eq_iInf, ramificationFiltration_Gv_eq_iInf]
  exact key

def isNaturalFiltration_of_isNaturalStage.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

/-! ## §4 ★★非空虚性 —— 仮説は空ではない(すべて無条件) -/

/-- ★★**抽象核(純群論)** —— 「内部自己同型でずれるだけ」の自己同型は
**正規部分群を保つ**。`map_eq_map_of_conj` の `Ψ = id` 版だが、`Subgroup.map_id`
との噛み合わせを避けて直に書く方が短い。 -/
theorem map_eq_self_of_conj {Γ : Type*} [Group Γ] (Φ : Γ ≃* Γ) (H : Subgroup Γ) [H.Normal]
    (h : ∀ g : Γ, ∃ c : Γ, Φ g = c * g * c⁻¹) : Subgroup.map (Φ : Γ →* Γ) H = H := by
  rw [map_eq_map_of_conj Φ (MulEquiv.refl Γ) H h]
  ext x
  simp

/-- ★★**`v ≤ 0` では段データの自然性は無条件に成り立つ**。

`absGalStage K N v = I_K ⊔ N`(`absGalStage_of_nonpos`)と
Corollary 1.3(`inertia_recoverable_real`、無条件)を合わせるだけ。

⇒ ★`IsNaturalStage` を `v > 0` に絞ったのは**弱め方ではない**。 -/
theorem map_absGalStage_of_nonpos {K K' : PAdicLocalField p}
    (Φ : ContinuousMulEquiv K.absGal K'.absGal) {N : Subgroup K.absGal}
    (hN : N ∈ openNormalBase K.absGal) {v : ℝ} (hv : v ≤ 0) :
    Subgroup.map ((Φ.toMulEquiv : K.absGal ≃* K'.absGal) : K.absGal →* K'.absGal)
        (absGalStage K N v)
      = absGalStage K' (Subgroup.map ((Φ.toMulEquiv : K.absGal ≃* K'.absGal) :
          K.absGal →* K'.absGal) N) v := by
  rw [absGalStage_of_nonpos K hN hv,
    absGalStage_of_nonpos K' (map_mem_openNormalBase Φ hN) hv, Subgroup.map_sup]
  congr 1
  rw [absInertia_eq_inertia, absInertia_eq_inertia]
  exact inertia_recoverable_real Φ

def map_absGalStage_of_nonpos.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Corollary 1.3", sectionId := "cor-1-3" }

/-- ★★`IsNaturalStage`(`v > 0` だけ)は**全実数版と同値**。
——`v ≤ 0` が無条件だから(`map_absGalStage_of_nonpos`)。 -/
theorem isNaturalStage_iff_forall :
    IsNaturalStage p ↔
      ∀ {K K' : PAdicLocalField p} (α : K.carrier ≃ₐ[ℚ_[p]] K'.carrier)
        ⦃N : Subgroup K.absGal⦄, N ∈ openNormalBase K.absGal → ∀ v : ℝ,
          Subgroup.map (((galContinuousMulEquiv α).toMulEquiv : K.absGal ≃* K'.absGal) :
              K.absGal →* K'.absGal) (absGalStage K N v)
            = absGalStage K' (Subgroup.map (((galContinuousMulEquiv α).toMulEquiv :
                K.absGal ≃* K'.absGal) : K.absGal →* K'.absGal) N) v := by
  constructor
  · intro h _ _ α N hN v
    rcases lt_or_ge 0 v with hv | hv
    · exact h α hN v hv
    · exact map_absGalStage_of_nonpos (galContinuousMulEquiv α) hN hv
  · intro h _ _ α N hN v _
    exact h α hN v

/-- ★★**内部自己同型では段データの自然性は無条件**——
`absGalStage K N v` が正規(`normal_absGalStage`)で `N` も正規だから。 -/
theorem map_absGalStage_inner (K : PAdicLocalField p) (c : K.absGal)
    {N : Subgroup K.absGal} (hN : N ∈ openNormalBase K.absGal) (v : ℝ) :
    Subgroup.map (((innerAbsGalEquiv K c).toMulEquiv : K.absGal ≃* K.absGal) :
        K.absGal →* K.absGal) (absGalStage K N v)
      = absGalStage K (Subgroup.map (((innerAbsGalEquiv K c).toMulEquiv : K.absGal ≃* K.absGal) :
          K.absGal →* K.absGal) N) v := by
  haveI := hN.2
  haveI := normal_absGalStage K hN v
  have hN' : Subgroup.map (((innerAbsGalEquiv K c).toMulEquiv : K.absGal ≃* K.absGal) :
      K.absGal →* K.absGal) N = N :=
    map_eq_self_of_conj _ _ (fun g => ⟨c, rfl⟩)
  rw [hN']
  exact map_eq_self_of_conj _ _ (fun g => ⟨c, rfl⟩)

/-- ★★**内部自己同型は `Γ_K^v` を保つ**(どの `RamificationFiltration` でも、無条件)。 -/
theorem map_gv_inner (RF : RamificationFiltration p) (K : PAdicLocalField p) (c : K.absGal)
    (v : ℝ) :
    Subgroup.map (((innerAbsGalEquiv K c).toMulEquiv : K.absGal ≃* K.absGal) :
        K.absGal →* K.absGal) (RF.Gv K v) = RF.Gv K v := by
  haveI := RF.isNormal K v
  exact map_eq_self_of_conj _ _ (fun g => ⟨c, rfl⟩)

/-- ★★★**`α = id` の場合の自然性は無条件に成り立つ**(どの `RamificationFiltration` でも)。

★これは自明ではない —— `galContinuousMulEquiv (AlgEquiv.refl)` は
`extendToClosure (AlgEquiv.refl)` の選択のせいで**恒等写像とは限らない**。
実際には `Γ_K` の内部自己同型であり(`galMulEquiv_conj_indep` を `ᾱ = RingEquiv.refl`
に当てる)、正規性で吸収される。 -/
theorem map_gv_refl (RF : RamificationFiltration p) (K : PAdicLocalField p) (v : ℝ) :
    Subgroup.map (((galContinuousMulEquiv (AlgEquiv.refl :
          K.carrier ≃ₐ[ℚ_[p]] K.carrier)).toMulEquiv : K.absGal ≃* K.absGal) :
        K.absGal →* K.absGal) (RF.Gv K v) = RF.Gv K v := by
  haveI := RF.isNormal K v
  have hfwd : ∀ x : K.carrier, (RingEquiv.refl K.closure) (algebraMap K.carrier K.closure x)
      = algebraMap K.carrier K.closure ((AlgEquiv.refl : K.carrier ≃ₐ[ℚ_[p]] K.carrier) x) := by
    intro x; rfl
  have hid : ∀ g : K.absGal,
      galMulEquivOf (AlgEquiv.refl : K.carrier ≃ₐ[ℚ_[p]] K.carrier)
        (RingEquiv.refl K.closure) hfwd g = g := by
    intro g; exact AlgEquiv.ext (fun _ => rfl)
  have hconj : ∀ g : K.absGal, ∃ c : K.absGal,
      (galContinuousMulEquiv (AlgEquiv.refl : K.carrier ≃ₐ[ℚ_[p]] K.carrier)).toMulEquiv g
        = c * g * c⁻¹ := by
    intro g
    obtain ⟨c, hc⟩ := galMulEquivOf_indep (AlgEquiv.refl : K.carrier ≃ₐ[ℚ_[p]] K.carrier)
      (extendToClosure (AlgEquiv.refl : K.carrier ≃ₐ[ℚ_[p]] K.carrier))
      (RingEquiv.refl K.closure)
      (extendToClosure_algebraMap (AlgEquiv.refl : K.carrier ≃ₐ[ℚ_[p]] K.carrier)) hfwd g
    rw [hid g] at hc
    exact ⟨c, hc⟩
  exact map_eq_self_of_conj _ _ hconj

def map_gv_refl.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Theorem 4.2", sectionId := "theorem-4-2" }

#print axioms map_eq_map_of_conj
#print axioms map_iInf_eq_of_bij
#print axioms map_mem_openNormalBase
#print axioms comap_mem_openNormalBase
#print axioms map_gv_eq_of_ext
#print axioms ramificationFiltration_Gv_eq_iInf
#print axioms isNaturalFiltration_of_isNaturalStage
#print axioms map_eq_self_of_conj
#print axioms map_absGalStage_of_nonpos
#print axioms isNaturalStage_iff_forall
#print axioms map_absGalStage_inner
#print axioms map_gv_inner
#print axioms map_gv_refl

end ABC3.Found.PGC
