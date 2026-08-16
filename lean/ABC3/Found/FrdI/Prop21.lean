import ABC3.Found.FrdI.Prop114

/-!
# [FrdI] Proposition 2.1 —— The Naive Frobenius Functor

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.44–p.45。

原文 (FrdI p.44):
> Proposition 2.1.
> (The Naive Frobenius Functor) Let d ∈N≥1. Then:

## ★この命題の規模(測定)

**3 条 (i)–(iii)**、主張の数は **7**:

| 条 | 主張の数 | 内容 |
|---|---|---|
| (i) | 3 | 関手であること / 同型を除いて一意 / 次数の合成則 `Ψ_{d₁} ∘ Ψ_{d₂} ≅ Ψ_{d₁d₂}` |
| (ii) | 2 | `𝒞 → 𝔽_Φ` に対する 1-両立性 / `𝒪^▷` 上で `d` 乗になること |
| (iii) | 2 | perfect 型 ⟺ `Ψ` が圏同値(両向き) |

## ★構成の要点

★**対象の側は選択が要る** —— `Definition 1.3, (ii)`(`frobDegSurj`)は
「次数 `d` の Frobenius 型射が**存在する**」としか言わないので、
`Ψ(A)` はその終域を**選んだ**もの。★**だから原文も「well-defined up to isomorphism」
と書いている。**

★**射の側は一意**である —— `Proposition 1.10, (i)`
(`prop_1_10_i_exists_given`)が `∃!` を与える。
★**関手性(`map_id` / `map_comp`)はこの一意性から出る**(可換図式を 2 通りに読むだけ)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★(i) —— naive Frobenius 関手

原文 (FrdI p.44):
> d2 ∈N≥1 on C is isomorphic to the naive Frobenius functor of degree d1 · d2

原文 (FrdI p.44):
> [well-defined up to isomorphism of functors] which we shall refer to as the naive
-/

section NaiveFrobenius

variable (F : FrobenioidCore P) (d : ℕ+)

/-- ★`Ψ(A)` —— 次数 `d` の Frobenius 型射の終域として**選んだ**対象。

★`Definition 1.3, (ii)` は存在しか言わないので、ここで選択が入る。 -/
noncomputable def nfObj (A : C) : C := (F.frobDegSurj A d).choose

/-- ★`α_A : A ⟶ Ψ(A)` —— 選んだ Frobenius 型射。 -/
noncomputable def nfHom (A : C) : A ⟶ nfObj P F d A :=
  (F.frobDegSurj A d).choose_spec.choose

theorem nfHom_frobType (A : C) : IsFrobeniusType P (nfHom P F d A) :=
  (F.frobDegSurj A d).choose_spec.choose_spec.1

theorem nfHom_degFr (A : C) : P.degFr (nfHom P F d A) = d :=
  (F.frobDegSurj A d).choose_spec.choose_spec.2

/-- ★`Ψ(φ)` —— `Proposition 1.10, (i)` の**一意な** `φ'`。 -/
noncomputable def nfMap {A B : C} (φ : A ⟶ B) : nfObj P F d A ⟶ nfObj P F d B :=
  (prop_1_10_i_exists_given P F φ (nfHom P F d A) (nfHom_frobType P F d A)
    (nfHom P F d B) (nfHom_frobType P F d B)
    (by simp [nfHom_degFr])).choose

/-- ★定義する四角形。 -/
theorem nfMap_sq {A B : C} (φ : A ⟶ B) :
    φ ≫ nfHom P F d B = nfHom P F d A ≫ nfMap P F d φ :=
  (prop_1_10_i_exists_given P F φ (nfHom P F d A) (nfHom_frobType P F d A)
    (nfHom P F d B) (nfHom_frobType P F d B)
    (by simp [nfHom_degFr])).choose_spec.1

/-- ★四角形を満たすものは `nfMap` に限る(`Proposition 1.10, (i)` の一意性)。 -/
theorem nfMap_uniq {A B : C} (φ : A ⟶ B) (ψ : nfObj P F d A ⟶ nfObj P F d B)
    (h : φ ≫ nfHom P F d B = nfHom P F d A ≫ ψ) : ψ = nfMap P F d φ :=
  prop_1_10_i_uniq P.totEpiC φ (nfHom P F d A) (nfHom P F d B) ψ (nfMap P F d φ)
    h (nfMap_sq P F d φ)

/-- ★★**[FrdI] Proposition 2.1, (i)** —— **naive Frobenius 関手** `Ψ : 𝒞 ⥤ 𝒞`。

★**関手性は `Proposition 1.10, (i)` の一意性だけから出る。** -/
noncomputable def naiveFrob : C ⥤ C where
  obj := nfObj P F d
  map := nfMap P F d
  map_id A := (nfMap_uniq P F d (𝟙 A) (𝟙 _) (by simp)).symm
  map_comp {A B E} f g := by
    refine (nfMap_uniq P F d (f ≫ g) (nfMap P F d f ≫ nfMap P F d g) ?_).symm
    have h1 := nfMap_sq P F d f
    have h2 := nfMap_sq P F d g
    conv_lhs => rw [Category.assoc, h2, ← Category.assoc, h1, Category.assoc]

@[simp] theorem naiveFrob_obj (A : C) : (naiveFrob P F d).obj A = nfObj P F d A := rfl

@[simp] theorem naiveFrob_map {A B : C} (φ : A ⟶ B) :
    (naiveFrob P F d).map φ = nfMap P F d φ := rfl

/-- ★★**`α` は自然変換 `𝟭 𝒞 ⟶ Ψ` をなす**。

★原文の四角形 `φ′ ◦ α = β ◦ φ` は**まさに自然性**である。 -/
noncomputable def naiveFrobUnit : 𝟭 C ⟶ naiveFrob P F d where
  app A := nfHom P F d A
  naturality {A B} φ := nfMap_sq P F d φ

end NaiveFrobenius

/-! ## ★(i) の後半 —— 選び方に依らないこと

原文 (FrdI p.44):
> [well-defined up to isomorphism of functors] which we shall refer to as the naive
-/

/-- ★★**選び方を変えても関手は同型** —— 「well-defined up to isomorphism of functors」。

★★**中身は `Definition 1.3, (ii)` の本質的一意性**(`frobDegUniq`) ——
2 つの選択 `α_A`・`α'_A` は同型 `θ_A` で繋がり、
★**その `θ` が自然変換になることは `Proposition 1.10, (i)` の一意性から出る。** -/
noncomputable def nfCompare (F F' : FrobenioidCore P) (d : ℕ+) (A : C) :
    nfObj P F d A ⟶ nfObj P F' d A :=
  (F.frobDegUniq A (nfObj P F d A) (nfObj P F' d A) (nfHom P F d A) (nfHom P F' d A)
    (nfHom_frobType P F d A) (nfHom_frobType P F' d A)
    (by simp [nfHom_degFr])).choose

theorem nfCompare_isIso (F F' : FrobenioidCore P) (d : ℕ+) (A : C) :
    IsIso (nfCompare P F F' d A) :=
  (F.frobDegUniq A (nfObj P F d A) (nfObj P F' d A) (nfHom P F d A) (nfHom P F' d A)
    (nfHom_frobType P F d A) (nfHom_frobType P F' d A)
    (by simp [nfHom_degFr])).choose_spec.1

theorem nfCompare_sq (F F' : FrobenioidCore P) (d : ℕ+) (A : C) :
    nfHom P F d A ≫ nfCompare P F F' d A = nfHom P F' d A :=
  (F.frobDegUniq A (nfObj P F d A) (nfObj P F' d A) (nfHom P F d A) (nfHom P F' d A)
    (nfHom_frobType P F d A) (nfHom_frobType P F' d A)
    (by simp [nfHom_degFr])).choose_spec.2

/-- ★`nfCompare` は自然 —— `Proposition 1.10, (i)` の一意性から。 -/
theorem nfCompare_naturality (F F' : FrobenioidCore P) (d : ℕ+) {A B : C} (φ : A ⟶ B) :
    nfMap P F d φ ≫ nfCompare P F F' d B = nfCompare P F F' d A ≫ nfMap P F' d φ := by
  refine prop_1_10_i_uniq P.totEpiC φ (nfHom P F d A) (nfHom P F' d B) _ _ ?_ ?_
  · rw [← Category.assoc, ← nfMap_sq P F d φ, Category.assoc, nfCompare_sq]
  · rw [← Category.assoc, nfCompare_sq, ← nfMap_sq P F' d φ]

noncomputable def naiveFrobIso (F F' : FrobenioidCore P) (d : ℕ+) :
    naiveFrob P F d ≅ naiveFrob P F' d :=
  NatIso.ofComponents
    (fun A => @asIso _ _ _ _ (nfCompare P F F' d A) (nfCompare_isIso P F F' d A))
    (fun φ => nfCompare_naturality P F F' d φ)

/-! ## ★(ii) —— `𝔽_Φ` 上の Frobenius 関手との 1-両立性

原文 (FrdI p.44):
> (ii) The functor Ψ of (i) is “1-compatible”, relative to C →FΦ, with the

原文 (FrdI p.44):
> Definition 1.1, (iii)] by the endomorphism of the functor Φ

原文 (FrdI p.44):
> by Ψ is given by raising to the d-th power.

★**`𝔽_Φ` 上の Frobenius 関手は `Φ` の自己準同型「`d` 倍」で決まる**
(`Definition 1.1, (iii)`)ので、1-両立性の中身は
★**`Ψ(φ)` の不変量が `(Base φ, d · α*(Div φ), degFr φ)` であること**である。
-/

/-- ★**(ii)** 次数は保たれる。 -/
theorem prop_2_1_ii_degFr (F : FrobenioidCore P) (d : ℕ+) {A B : C} (φ : A ⟶ B) :
    P.degFr (nfMap P F d φ) = P.degFr φ :=
  (prop_1_10_i_degFr_phi_eq P (by simp [nfHom_degFr]) (nfMap_sq P F d φ)).symm

/-- ★★**(ii)** `Div` は「`d` 倍」で移る —— ★**これが `𝔽_Φ` 上の
Frobenius 関手(＝`Φ` の `d` 倍)との 1-両立性の中身**である。 -/
theorem prop_2_1_ii_Div (F : FrobenioidCore P) (d : ℕ+) {A B : C} (φ : A ⟶ B) :
    haveI : IsIso (P.Base (nfHom P F d A)) := (nfHom_frobType P F d A).2
    P.Div (nfMap P F d φ)
      = ((d : ℕ+) : ℕ) • Φ.map (inv (P.Base (nfHom P F d A))) (P.Div φ) := by
  haveI : IsIso (P.Base (nfHom P F d A)) := (nfHom_frobType P F d A).2
  have h := prop_1_10_i_Div_formula' P φ (nfHom P F d A) (nfHom P F d B) (nfMap P F d φ)
    (nfHom_frobType P F d A) (nfHom_frobType P F d B).1.2 (nfMap_sq P F d φ)
  rw [h, nfHom_degFr]

include P in
/-- ★★★**(ii) の後半** —— `A` が Frobenius-normalized で `α` が底恒等なら、
`Ψ` が誘導する `𝒪^▷(A) → 𝒪^▷(A)` は **`d` 乗**である。

★**中身は `Definition 1.2, (iv)` の Frobenius-normalized の等式そのもの** ——
`ζ ≫ α^d = α ≫ ζ` と、`Proposition 1.10, (i)` の四角形 `α ≫ ζ = ζ ≫ α'` を
突き合わせ、★**`𝒞` が totally epimorphic なので `ζ` を左から消約する。** -/
theorem prop_2_1_ii_pow {A : C} (d : ℕ+) (hfn : IsFrobeniusNormalized P A)
    (ζ : End A) (hζb : IsBaseIdentity P ζ) (hζd : P.degFr (ζ : A ⟶ A) = d)
    (α α' : End A) (hα : α ∈ OTri P A)
    (hsq : (α : A ⟶ A) ≫ (ζ : A ⟶ A) = (ζ : A ⟶ A) ≫ (α' : A ⟶ A)) :
    α' = α ^ ((d : ℕ+) : ℕ) := by
  have h := hfn ζ hζb α hα
  rw [hζd] at h
  haveI : Epi (ζ : A ⟶ A) := P.totEpiC _ _ _
  refine (cancel_epi (ζ : A ⟶ A)).mp ?_
  show (ζ : A ⟶ A) ≫ (α' : A ⟶ A) = (ζ : A ⟶ A) ≫ ((α ^ ((d : ℕ+) : ℕ) : End A) : A ⟶ A)
  rw [← hsq, h]

/-! ## ★(iii) —— perfect 型 ⟺ `Ψ` が圏同値

原文 (FrdI p.44):
> (iii) C is of perfect type if and only if Ψ is an equivalence of categories.

★★**我々の読み(量化子について)**: 原文は命題の冒頭で `d ∈ ℕ≥1` を固定するが、
`perfect object` の定義(`Definition 1.2, (iv)`)は **すべての `n ∈ ℕ≥1`** を量化する
(「for every n ∈ N≥1, ... every B ∈ Ob(C) base-isomorphic to C appears as the codomain
of a morphism of Frobenius type of Frobenius degree n」)。
★したがって `⟸` は **1 つの `d`** だけからは出ない —— `Ψ_d` が圏同値でも
`n ≠ d` についての perfect 性は言えない。
★**そこで `⟸` は「すべての `d` について `Ψ_d` が圏同値」の形で述べる。**
`⟹` は `d` を固定したままで成り立つ(`prop_2_1_iii_mp`)。

★**原文の証明の筋をそのまま辿る**:
- **本質的全射性**: perfect の第 1 条件そのもの
- **忠実性・充満性**: `Definition 1.3, (iv), (a)` の分解と `𝒞` の totally epimorphic 性で
  **linear の場合**に帰着し、`Definition 1.3, (i), (c)`(pull-back)で
  さらに **pre-step の場合**に帰着し、そこは perfect の第 2 条件(`∃!`)が与える
-/

section PerfectType

variable (F : FrobenioidCore P)

/-! ### ★道具 —— pull-back の普遍性を 2 つに割る -/

include P in
/-- ★**pull-back に沿った消約** —— 合成と底が一致すれば射は一致する。 -/
theorem isPullBack_inj {Z B : C} {α : Z ⟶ B} (h : IsPullBack P α) {X : C} {f g : X ⟶ Z}
    (h1 : f ≫ α = g ≫ α) (h2 : P.Base f = P.Base g) : f = g :=
  (h X).1 (Subtype.ext (Prod.ext h1 h2))

include P in
/-- ★**pull-back に沿った持ち上げ** —— 底を指定して一意に持ち上がる。 -/
theorem isPullBack_lift {Z B : C} {α : Z ⟶ B} (h : IsPullBack P α) {X : C} (g : X ⟶ B)
    (b : (P.toElem.obj X).base ⟶ (P.toElem.obj Z).base)
    (hb : P.Base g = b ≫ P.Base α) :
    ∃ f : X ⟶ Z, f ≫ α = g ∧ P.Base f = b := by
  obtain ⟨f, hf⟩ := (h X).2 ⟨(g, b), hb⟩
  have hfv := congrArg Subtype.val hf
  exact ⟨f, congrArg Prod.fst hfv, congrArg Prod.snd hfv⟩

include P F in
/-- ★**linear な射は「pre-step ≫ pull-back」に分解する**。

★`Definition 1.3, (iv), (a)` の 3 分解で Frobenius 部分の次数が 1 になり、
`Proposition 1.4, (iii)` によりそれは**同型**だから、pre-step に吸収できる。 -/
theorem linear_factor {A B : C} (v : A ⟶ B) (hv : IsLinear P v) :
    ∃ (Z : C) (β : A ⟶ Z) (α : Z ⟶ B), v = β ≫ α ∧ IsPreStep P β ∧ IsPullBack P α := by
  obtain ⟨X, Y, γ, β, α, hv', hγ, hβ, hα⟩ := F.arbFactor v
  have hdγ : P.degFr γ = 1 := by
    have h1 : ((P.degFr α : ℕ+) : ℕ)
        * (((P.degFr β : ℕ+) : ℕ) * ((P.degFr γ : ℕ+) : ℕ)) = 1 := by
      rw [← PNat.mul_coe, ← PNat.mul_coe, ← P.degFr_comp, ← P.degFr_comp,
        Category.assoc, ← hv', show P.degFr v = 1 from hv]
      rfl
    have hdvd : ((P.degFr γ : ℕ+) : ℕ) ∣ 1 :=
      ⟨((P.degFr α : ℕ+) : ℕ) * ((P.degFr β : ℕ+) : ℕ), by rw [← h1]; ring⟩
    exact PNat.coe_injective (by rw [Nat.dvd_one.mp hdvd]; rfl)
  haveI : IsIso γ := prop_1_4_iii P F γ hγ.1 ⟨hdγ, hγ.2⟩
  exact ⟨Y, γ ≫ β, α, by rw [hv', Category.assoc],
    IsPreStep.comp P (isPreStep_of_isIso P γ) hβ, hα⟩

/-! ### ★道具 —— `Ψ` が保つもの・反映するもの -/

variable (d : ℕ+)

/-- ★`Ψ(φ)` の底は `Base φ` を両端の同型で挟んだもの。 -/
theorem nfMap_base {A B : C} (φ : A ⟶ B) :
    haveI : IsIso (P.Base (nfHom P F d A)) := (nfHom_frobType P F d A).2
    P.Base (nfMap P F d φ)
      = inv (P.Base (nfHom P F d A)) ≫ P.Base φ ≫ P.Base (nfHom P F d B) := by
  haveI : IsIso (P.Base (nfHom P F d A)) := (nfHom_frobType P F d A).2
  have h := congrArg P.Base (nfMap_sq P F d φ)
  rw [P.Base_comp, P.Base_comp] at h
  rw [h, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]

/-- ★`Ψ` は合成を保つ(関手性の再掲、`rw` 用)。 -/
theorem nfMap_comp {A B E : C} (f : A ⟶ B) (g : B ⟶ E) :
    nfMap P F d (f ≫ g) = nfMap P F d f ≫ nfMap P F d g :=
  (naiveFrob P F d).map_comp f g

/-- ★`Ψ` は底について忠実 —— `Base (Ψ φ)` が決まれば `Base φ` が決まる。 -/
theorem nfMap_base_inj {A B : C} {φ φ' : A ⟶ B}
    (h : P.Base (nfMap P F d φ) = P.Base (nfMap P F d φ')) : P.Base φ = P.Base φ' := by
  haveI : IsIso (P.Base (nfHom P F d B)) := (nfHom_frobType P F d B).2
  have h1 := congrArg P.Base (nfMap_sq P F d φ)
  have h2 := congrArg P.Base (nfMap_sq P F d φ')
  rw [P.Base_comp, P.Base_comp] at h1
  rw [P.Base_comp, P.Base_comp] at h2
  rw [h] at h1
  exact (cancel_mono (P.Base (nfHom P F d B))).mp (h1.trans h2.symm)

theorem nfMap_frobType {A B : C} (φ : A ⟶ B) (h : IsFrobeniusType P φ) :
    IsFrobeniusType P (nfMap P F d φ) :=
  prop_1_10_i_frobType_of P F (nfHom_frobType P F d A) (nfHom_frobType P F d B)
    (nfMap_sq P F d φ) h

theorem nfMap_preStep {A B : C} (φ : A ⟶ B) (h : IsPreStep P φ) :
    IsPreStep P (nfMap P F d φ) :=
  prop_1_10_i_preStep_of P (nfHom_frobType P F d A) (nfHom_frobType P F d B)
    (by simp [nfHom_degFr]) (nfMap_sq P F d φ) h

theorem nfMap_pullBack {A B : C} (φ : A ⟶ B) (h : IsPullBack P φ) :
    IsPullBack P (nfMap P F d φ) :=
  prop_1_10_i_pullBack_of P F (nfHom_frobType P F d A) (nfHom_frobType P F d B)
    (by simp [nfHom_degFr]) (nfMap_sq P F d φ) h

/-- ★`Ψ` は pre-step を**反映**する。 -/
theorem nfMap_preStep_of {A B : C} (φ : A ⟶ B) (h : IsPreStep P (nfMap P F d φ)) :
    IsPreStep P φ := by
  haveI : IsIso (P.Base (nfHom P F d A)) := (nfHom_frobType P F d A).2
  haveI : IsIso (P.Base (nfHom P F d B)) := (nfHom_frobType P F d B).2
  haveI : IsIso (P.Base (nfMap P F d φ)) := h.2
  refine ⟨?_, ?_⟩
  · show P.degFr φ = 1
    rw [← prop_2_1_ii_degFr P F d φ]; exact h.1
  · show IsIso (P.Base φ)
    have h1 := congrArg P.Base (nfMap_sq P F d φ)
    rw [P.Base_comp, P.Base_comp] at h1
    have h2 : P.Base φ
        = (P.Base (nfHom P F d A) ≫ P.Base (nfMap P F d φ))
          ≫ inv (P.Base (nfHom P F d B)) := by
      rw [← h1, Category.assoc, IsIso.hom_inv_id, Category.comp_id]
    rw [h2]; infer_instance

/-- ★`A` と `Ψ(A)` は base-isomorphic。 -/
theorem baseIsomorphic_nfObj (A : C) : BaseIsomorphic P A (nfObj P F d A) := by
  haveI : IsIso (P.Base (nfHom P F d A)) := (nfHom_frobType P F d A).2
  exact ⟨asIso (P.Base (nfHom P F d A))⟩

/-! ### ★`⟹` —— perfect 型なら `Ψ_d` は圏同値 -/

/-- ★**本質的全射性** —— perfect の第 1 条件そのもの。 -/
theorem naiveFrob_essSurj (hpt : IsOfPerfectType P) : (naiveFrob P F d).EssSurj := by
  constructor
  intro B
  obtain ⟨B₀, φ, hφF, hφd⟩ := (hpt B d).1 B ⟨Iso.refl _⟩
  obtain ⟨β, hβiso, -⟩ := F.frobDegUniq B₀ (nfObj P F d B₀) B
    (nfHom P F d B₀) φ (nfHom_frobType P F d B₀) hφF (by rw [nfHom_degFr, hφd])
  haveI : IsIso β := hβiso
  exact ⟨B₀, ⟨@asIso _ _ _ _ β hβiso⟩⟩

/-- ★★**忠実性** —— 3 分解で linear に落とし、pull-back で pre-step に落として、
perfect の第 2 条件の**一意性**で締める。 -/
theorem naiveFrob_faithful (hpt : IsOfPerfectType P) : (naiveFrob P F d).Faithful := by
  constructor
  intro A B φ φ' hmap
  show φ = φ'
  have hmap' : nfMap P F d φ = nfMap P F d φ' := hmap
  -- 手 1: Frobenius 部分を揃えて linear に落とす
  obtain ⟨X, Y, γ, β, α, hφ, hγ, hβ, hα⟩ := F.arbFactor φ
  obtain ⟨X', Y', γ', β', α', hφ', hγ', hβ', hα'⟩ := F.arbFactor φ'
  have hlin : IsLinear P (β ≫ α) :=
    IsLinear.comp P hβ.1 (F.pullBackLB α hα).2
  have hlin' : IsLinear P (β' ≫ α') :=
    IsLinear.comp P hβ'.1 (F.pullBackLB α' hα').2
  have hdd : P.degFr γ = P.degFr γ' := by
    have e1 : P.degFr φ = P.degFr γ := by
      rw [hφ, P.degFr_comp, show P.degFr (β ≫ α) = 1 from hlin, one_mul]
    have e2 : P.degFr φ' = P.degFr γ' := by
      rw [hφ', P.degFr_comp, show P.degFr (β' ≫ α') = 1 from hlin', one_mul]
    rw [← e1, ← e2, ← prop_2_1_ii_degFr P F d φ, ← prop_2_1_ii_degFr P F d φ', hmap']
  obtain ⟨θ, hθ, hθe⟩ := F.frobDegUniq A X X' γ γ' hγ hγ' hdd
  haveI := hθ
  set v : X ⟶ B := β ≫ α with hvdef
  set w : X ⟶ B := θ ≫ β' ≫ α' with hwdef
  have hφw : φ' = γ ≫ w := by rw [hwdef, ← Category.assoc, hθe, hφ']
  have hwlin : IsLinear P w := IsLinear.comp P (isLinear_of_isIso P θ) hlin'
  -- 手 2: `Ψ(γ)` は epi なので linear 部分に落ちる
  have hvw : nfMap P F d v = nfMap P F d w := by
    haveI : Epi (nfMap P F d γ) := P.totEpiC _ _ _
    refine (cancel_epi (nfMap P F d γ)).mp ?_
    have e1 : nfMap P F d γ ≫ nfMap P F d v = nfMap P F d φ := by
      rw [← nfMap_comp, ← hφ]
    have e2 : nfMap P F d γ ≫ nfMap P F d w = nfMap P F d φ' := by
      rw [← nfMap_comp, ← hφw]
    rw [e1, e2, hmap']
  -- 手 3: pull-back で pre-step に落とす
  obtain ⟨Z, β₁, α₁, hv1, hβ₁, hα₁⟩ := linear_factor P F v hlin
  have hbvw : P.Base v = P.Base w := nfMap_base_inj P F d (by rw [hvw])
  obtain ⟨u, hu1, hu2⟩ := isPullBack_lift P hα₁ w (P.Base β₁) (by
    rw [← hbvw, hv1, P.Base_comp])
  have hulin : IsLinear P u := by
    have : P.degFr w = P.degFr α₁ * P.degFr u := by rw [← hu1, P.degFr_comp]
    rw [show P.degFr w = 1 from hwlin,
      show P.degFr α₁ = 1 from (F.pullBackLB α₁ hα₁).2, one_mul] at this
    exact this.symm
  have hus : IsPreStep P u := ⟨hulin, by show IsIso (P.Base u); rw [hu2]; exact hβ₁.2⟩
  -- 手 4: `Ψ(α₁)` の pull-back 性で pre-step 部分を分離する
  have hβu : nfMap P F d β₁ = nfMap P F d u := by
    refine isPullBack_inj P (nfMap_pullBack P F d α₁ hα₁) ?_ ?_
    · have e1 : nfMap P F d β₁ ≫ nfMap P F d α₁ = nfMap P F d v := by
        rw [← nfMap_comp, ← hv1]
      have e2 : nfMap P F d u ≫ nfMap P F d α₁ = nfMap P F d w := by
        rw [← nfMap_comp, hu1]
      rw [e1, e2, hvw]
    · rw [nfMap_base, nfMap_base, hu2]
  -- 手 5: perfect の一意性
  have hbi : BaseIsomorphic P Z X := by
    haveI : IsIso (P.Base β₁) := hβ₁.2
    exact ⟨(asIso (P.Base β₁)).symm⟩
  obtain ⟨-, huniq⟩ := hpt X d
  have hβ₁u : β₁ = u := by
    obtain ⟨s, -, hs⟩ := huniq X (nfObj P F d X) Z (nfObj P F d Z)
      (nfHom P F d X) (nfHom P F d Z) (nfHom_frobType P F d X) (nfHom_degFr P F d X)
      (nfHom_frobType P F d Z) (nfHom_degFr P F d Z) ⟨Iso.refl _⟩ hbi
      (nfMap P F d β₁) (nfMap_preStep P F d β₁ hβ₁)
    rw [hs β₁ ⟨hβ₁, (nfMap_sq P F d β₁).symm⟩,
      hs u ⟨hus, by rw [hβu]; exact (nfMap_sq P F d u).symm⟩]
  rw [hφ, hφw, hv1, hβ₁u, hu1]

/-- ★★**充満性(linear の場合)** —— `Ψ(A) ⟶ Ψ(B)` の linear な射は
`A ⟶ B` の linear な射から来る。

★`Definition 1.3, (i), (c)`(`plBk_realize`)で pull-back 部分を `𝒞` へ降ろし、
残る pre-step 部分は perfect の第 2 条件の**存在**が与える。 -/
theorem naiveFrob_full_linear (hpt : IsOfPerfectType P) {A B : C}
    (v : nfObj P F d A ⟶ nfObj P F d B) (hvlin : IsLinear P v) :
    ∃ u : A ⟶ B, IsLinear P u ∧ nfMap P F d u = v := by
  haveI hbB : IsIso (P.Base (nfHom P F d B)) := (nfHom_frobType P F d B).2
  -- 手 1: `v = β₁ ≫ α₁`
  obtain ⟨Z₁, β₁, α₁, hv1, hβ₁, hα₁⟩ := linear_factor P F v hvlin
  -- 手 2: `α₁` の底を `B` の上に降ろして pull-back `α₀` を作る
  obtain ⟨Z₀, α₀, hα₀, θ₀, hθ₀⟩ :=
    plBk_realize P F B (P.Base α₁ ≫ inv (P.Base (nfHom P F d B)))
  haveI hbZ₀ : IsIso (P.Base (nfHom P F d Z₀)) := (nfHom_frobType P F d Z₀).2
  have hα₀' : IsPullBack P (nfMap P F d α₀) := nfMap_pullBack P F d α₀ hα₀
  -- `Ψ(α₀)` の底は `b ≫ Base α₁`(`b` は同型)
  have hbase : P.Base (nfMap P F d α₀)
      = (inv (P.Base (nfHom P F d Z₀)) ≫ θ₀.hom) ≫ P.Base α₁ := by
    rw [nfMap_base, hθ₀, Category.assoc, Category.assoc, Category.assoc,
      IsIso.inv_hom_id, Category.comp_id]
  -- 手 3: 両側の普遍性で `Ψ(Z₀) ≅ Z₁` を作る
  obtain ⟨ρ, hρ1, hρ2⟩ := isPullBack_lift P hα₁ (nfMap P F d α₀)
    (inv (P.Base (nfHom P F d Z₀)) ≫ θ₀.hom) hbase
  obtain ⟨σ, hσ1, hσ2⟩ := isPullBack_lift P hα₀' α₁
    (inv (inv (P.Base (nfHom P F d Z₀)) ≫ θ₀.hom)) (by
      rw [hbase, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp])
  haveI hρiso : IsIso ρ := by
    refine ⟨σ, ?_, ?_⟩
    · refine isPullBack_inj P hα₀' ?_ ?_
      · rw [Category.assoc, hσ1, hρ1, Category.id_comp]
      · rw [P.Base_comp, hρ2, hσ2, IsIso.hom_inv_id, P.Base_id]
    · refine isPullBack_inj P hα₁ ?_ ?_
      · rw [Category.assoc, hρ1, hσ1, Category.id_comp]
      · rw [P.Base_comp, hσ2, hρ2, IsIso.inv_hom_id, P.Base_id]
  -- 手 4: `v = (β₁ ≫ ρ⁻¹) ≫ Ψ(α₀)`
  have hβ₂ : IsPreStep P (β₁ ≫ inv ρ) :=
    IsPreStep.comp P hβ₁ (isPreStep_of_isIso P (inv ρ))
  have hv2 : v = (β₁ ≫ inv ρ) ≫ nfMap P F d α₀ := by
    rw [← hρ1, ← Category.assoc, Category.assoc β₁, IsIso.inv_hom_id, Category.comp_id]
    exact hv1
  -- 手 5: perfect の第 2 条件で pre-step 部分を降ろす
  have hbi : BaseIsomorphic P Z₀ A := by
    haveI hbA : IsIso (P.Base (nfHom P F d A)) := (nfHom_frobType P F d A).2
    haveI hbβ₂ : IsIso (P.Base (β₁ ≫ inv ρ)) := hβ₂.2
    exact ⟨((asIso (P.Base (nfHom P F d Z₀))).trans
      (asIso (P.Base (β₁ ≫ inv ρ))).symm).trans
      (asIso (P.Base (nfHom P F d A))).symm⟩
  obtain ⟨-, hex⟩ := hpt A d
  obtain ⟨s, ⟨hs1, hs2⟩, -⟩ := hex A (nfObj P F d A) Z₀ (nfObj P F d Z₀)
    (nfHom P F d A) (nfHom P F d Z₀) (nfHom_frobType P F d A) (nfHom_degFr P F d A)
    (nfHom_frobType P F d Z₀) (nfHom_degFr P F d Z₀) ⟨Iso.refl _⟩ hbi
    (β₁ ≫ inv ρ) hβ₂
  refine ⟨s ≫ α₀, IsLinear.comp P hs1.1 (F.pullBackLB α₀ hα₀).2, ?_⟩
  rw [nfMap_comp, ← nfMap_uniq P F d s (β₁ ≫ inv ρ) hs2.symm, ← hv2]

/-- ★★**充満性** —— Frobenius 部分を `Definition 1.3, (ii)` で揃えてから
linear の場合(`naiveFrob_full_linear`)に帰着する。 -/
theorem naiveFrob_full (hpt : IsOfPerfectType P) : (naiveFrob P F d).Full := by
  constructor
  intro A B g
  -- 手 1: `g = γ ≫ (linear)`
  obtain ⟨W, Y, γ, β, α, hg, hγ, hβ, hα⟩ := F.arbFactor g
  obtain ⟨A₁, γA, hγA, hγAd⟩ := F.frobDegSurj A (P.degFr γ)
  obtain ⟨θ, hθ, hθe⟩ := F.frobDegUniq (nfObj P F d A) W (nfObj P F d A₁)
    γ (nfMap P F d γA) hγ (nfMap_frobType P F d γA hγA)
    (by rw [prop_2_1_ii_degFr]; exact hγAd.symm)
  haveI := hθ
  -- `g = Ψ(γA) ≫ v`
  have hgv : g = nfMap P F d γA ≫ (inv θ ≫ β ≫ α) := by
    rw [← hθe, Category.assoc, ← Category.assoc θ, IsIso.hom_inv_id, Category.id_comp]
    exact hg
  have hvlin : IsLinear P (inv θ ≫ β ≫ α) :=
    IsLinear.comp P (isLinear_of_isIso P (inv θ))
      (IsLinear.comp P hβ.1 (F.pullBackLB α hα).2)
  -- 手 2: linear の場合へ
  obtain ⟨u, -, hu⟩ := naiveFrob_full_linear P F d hpt _ hvlin
  exact ⟨γA ≫ u, by rw [naiveFrob_map, nfMap_comp, hu]; exact hgv.symm⟩

/-- ★★★**[FrdI] Proposition 2.1, (iii) の `⟹`** —— perfect 型なら `Ψ_d` は圏同値。 -/
theorem prop_2_1_iii_mp (hpt : IsOfPerfectType P) : (naiveFrob P F d).IsEquivalence := by
  haveI := naiveFrob_full P F d hpt
  haveI := naiveFrob_faithful P F d hpt
  haveI := naiveFrob_essSurj P F d hpt
  exact ⟨inferInstance, inferInstance, inferInstance⟩

/-- ★★★**[FrdI] Proposition 2.1, (iii) の `⟸`** —— **すべての次数**で `Ψ` が
圏同値なら perfect 型。

★★**量化子について**(冒頭の注記): `perfect object`(`Definition 1.2, (iv)`)は
**すべての `n ∈ ℕ≥1`** を量化するので、1 つの `d` だけでは足りない。 -/
theorem prop_2_1_iii_mpr (h : ∀ n : ℕ+, (naiveFrob P F n).IsEquivalence) :
    IsOfPerfectType P := by
  intro A n
  haveI := (h n).full
  haveI := (h n).faithful
  haveI := (h n).essSurj
  refine ⟨?_, ?_⟩
  · -- ★第 1 条件 —— 本質的全射性そのもの
    intro B _
    obtain ⟨B₀, ⟨e⟩⟩ := Functor.EssSurj.mem_essImage (F := naiveFrob P F n) B
    have e' : nfObj P F n B₀ ≅ B := e
    refine ⟨B₀, nfHom P F n B₀ ≫ e'.hom,
      IsFrobeniusType.comp P F (nfHom_frobType P F n B₀)
        (isFrobeniusType_of_isIso P e'.hom), ?_⟩
    rw [P.degFr_comp, show P.degFr e'.hom = 1 from isLinear_of_isIso P e'.hom,
      one_mul, nfHom_degFr]
  · -- ★第 2 条件 —— 充満性(存在)と忠実性(一意性)
    intro B₁ B₁' B₂ B₂' φ₁ φ₂ hφ₁ hd₁ hφ₂ hd₂ _ _ ψ' hψ'
    obtain ⟨θ₁, hθ₁, hθ₁e⟩ := F.frobDegUniq B₁ B₁' (nfObj P F n B₁) φ₁ (nfHom P F n B₁)
      hφ₁ (nfHom_frobType P F n B₁) (by rw [hd₁, nfHom_degFr])
    obtain ⟨θ₂, hθ₂, hθ₂e⟩ := F.frobDegUniq B₂ B₂' (nfObj P F n B₂) φ₂ (nfHom P F n B₂)
      hφ₂ (nfHom_frobType P F n B₂) (by rw [hd₂, nfHom_degFr])
    haveI hθ₁' : IsIso θ₁ := hθ₁
    haveI hθ₂' : IsIso θ₂ := hθ₂
    have hψ'' : IsPreStep P (inv θ₁ ≫ ψ' ≫ θ₂) :=
      IsPreStep.comp P (isPreStep_of_isIso P (inv θ₁))
        (IsPreStep.comp P hψ' (isPreStep_of_isIso P θ₂))
    obtain ⟨ψ, hψ⟩ :=
      (naiveFrob P F n).map_surjective (X := B₁) (Y := B₂)
        (show nfObj P F n B₁ ⟶ nfObj P F n B₂ from inv θ₁ ≫ ψ' ≫ θ₂)
    have hψ0 : nfMap P F n ψ = inv θ₁ ≫ ψ' ≫ θ₂ := hψ
    -- ★四角形を `φ₁`・`φ₂` の言葉に戻す
    have hback : ∀ (ξ : B₁ ⟶ B₂), (φ₁ ≫ ψ' = ξ ≫ φ₂)
        ↔ (ξ ≫ nfHom P F n B₂ = nfHom P F n B₁ ≫ (inv θ₁ ≫ ψ' ≫ θ₂)) := by
      intro ξ
      constructor
      · intro hs
        rw [← hθ₁e, ← hθ₂e, ← Category.assoc, ← hs]
        simp
      · intro hs
        rw [← hθ₁e, ← hθ₂e] at hs
        have h2 : (ξ ≫ φ₂) ≫ θ₂ = (φ₁ ≫ ψ') ≫ θ₂ := by
          rw [Category.assoc, Category.assoc, hs]; simp
        exact ((cancel_mono θ₂).mp h2).symm
    have hsq : φ₁ ≫ ψ' = ψ ≫ φ₂ :=
      (hback ψ).mpr (by rw [← hψ0]; exact nfMap_sq P F n ψ)
    refine ⟨ψ, ⟨nfMap_preStep_of P F n ψ (by rw [hψ0]; exact hψ''), hsq⟩, ?_⟩
    rintro ψ₂ ⟨-, hψ₂sq⟩
    refine (naiveFrob P F n).map_injective ?_
    show nfMap P F n ψ₂ = nfMap P F n ψ
    rw [hψ0, nfMap_uniq P F n ψ₂ (inv θ₁ ≫ ψ' ≫ θ₂) ((hback ψ₂).mp hψ₂sq)]

/-- ★★★**[FrdI] Proposition 2.1, (iii)**。

原文 (FrdI p.44):
> (iii) C is of perfect type if and only if Ψ is an equivalence of categories.

★★**量化子は原文どおりに読めない**(冒頭の注記): `Ψ` は次数 `d` ごとに定まるが、
`perfect` は**すべての次数**についての条件なので、`⟸` は「すべての `d`」を要する。
★`⟹` の側は `d` を固定したままで成り立つ(`prop_2_1_iii_mp`)。 -/
theorem prop_2_1_iii : IsOfPerfectType P ↔ ∀ n : ℕ+, (naiveFrob P F n).IsEquivalence :=
  ⟨fun hpt n => prop_2_1_iii_mp P F n hpt, prop_2_1_iii_mpr P F⟩

end PerfectType

end ABC3.Found.FrdI
