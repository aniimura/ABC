import ABC3.Found.FrdI.Prop22

/-!
# [FrdI] Definition 2.3 —— Characteristic Splitting

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.47。

原文 (FrdI p.47):
> We shall refer to as a characteristic splitting on C a subfunctor

原文 (FrdI p.47):
> following properties hold: (a) for every A

## ★この定義の中身(測定)

**主張は 4 つ**:

| # | 内容 |
|---|---|
| 1 | `τ` は `𝒪^▷(−) : (𝒞^istr)^lin → Mon` の**部分関手**である |
| 2 | (a) 各 `A` で `τ(A)` が `𝒪^▷(A)^char` へ**全単射**に写る |
| 3 | (a) ゆえに **`𝒪^×(A) × τ(A) ≅ 𝒪^▷(A)`** の分裂が `A` に関手的に定まる |
| 4 | (b) isotropic hull `A → A^istr` について `τ(A^istr)` が `𝒪^▷(A)` の像に入る |

★**3 は 2 から導く** —— 原文が「hence determines」と書くとおり。
★導出には `Proposition 2.2, (iii)` の**単射性**(`otri_div_eq_iff`)が要る:
`Div` が等しい 2 元は**単元 1 つ分しか違わない**。

## ★`𝒪^▷(A)^char` の扱い

`𝒪^▷(A)` は一般に**非可換**なので、`§0` の商モノイド `M/M^±` を
そのまま作るのは重い。★**代わりに `Div` の言葉で書く** ——
`Proposition 2.2, (iii)`(`otri_isUnit_iff_Div_zero` と `otri_div_eq_iff`)により

* `Div x = 0 ⟺ x ∈ 𝒪^×(A)`
* `Div x = Div y ⟺ x と y は単元 1 つ分しか違わない`

なので、★**`Div` の値がちょうど `𝒪^▷(A)^char` の元を指定する。**
「`τ(A)` が `𝒪^▷(A)^char` へ全単射に写る」は
「各 `x ∈ 𝒪^▷(A)` に対し `Div t = Div x` なる `t ∈ τ(A)` が**ただ 1 つ**」と同値である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P)

/-- ★★★**[FrdI] Definition 2.3** —— `𝒞` 上の **characteristic splitting**。

原文 (FrdI p.47):
> of Proposition 2.2, (ii), such that the

★`τ` は各対象に `End A` の部分モノイドを与える族として持つ。
★**条件 1・2 が `(𝒞^istr)^lin` 上の部分関手であること**、
★**条件 3 が (a)**、★**条件 4 が (b)** である。 -/
structure IsCharacteristicSplitting (τ : ∀ A : C, Submonoid (End A)) : Prop where
  /-- **部分関手 1** —— `τ(A) ⊆ 𝒪^▷(A)`。 -/
  le_otri : ∀ A : C, τ A ≤ OTri P A
  /-- **部分関手 2** —— `(𝒞^istr)^lin` の射(＝ isotropic 対象からの linear 射)で保たれる。 -/
  map_mem : ∀ {A B : C} (hA : IsIsotropic P A) {φ : A ⟶ B} (hl : IsLinear P φ)
    (t : OTri P B), ((t : End B)) ∈ τ B → ((otriLin P F hA hl t : End A)) ∈ τ A
  /-- **(a)** `τ(A)` は `𝒪^▷(A)^char` へ**全単射**に写る。

  原文 (FrdI p.47):
  > , hence determines a splitting of monoids
  -/
  charBij : ∀ A : C, IsIsotropic P A → ∀ x : OTri P A,
    ∃! t : τ A, P.Div (((t : End A)) : A ⟶ A) = P.Div (((x : End A)) : A ⟶ A)
  /-- **(b)** isotropic hull `A → A^istr` について、`τ(A^istr)` は
  `Proposition 2.2, (iv)` の単射による `𝒪^▷(A)` の像に入る。

  原文 (FrdI p.47):
  > via the natural injection of Proposition 2.2,
  -/
  hullMem : ∀ {A B : C} {φ : A ⟶ B} (hφ : IsIsotropicHull P φ) (t : End B), t ∈ τ B →
    ∃ s : OTri P A, hullOTriHom P φ hφ ((s : End A)) = t
  /-- **(b) の言い換え —— 非 isotropic な `A` での `τ(A)` は引き戻しである。**

  ★★**これは条件の追加ではない。** 原文の `τ` は
  `τ : (𝒞^istr)^lin → Mon` と、**`𝒞^istr` の上でのみ**定義された部分関手である
  (`Definition 2.3` 冒頭)。★したがって非 isotropic な `A` での `τ(A)` は
  原文では**この引き戻しが定義**であり、我々が `τ` を全対象で与えている以上、
  その対応を条件として書かねばならない。

  ★`A` が isotropic なら hull は同型なので、条件は自動的に満たされる。

  原文 (FrdI p.47):
  > via the natural injection of Proposition 2.2,
  -/
  hullPullback : ∀ {A B : C} {φ : A ⟶ B} (hφ : IsIsotropicHull P φ) (s : End A),
    s ∈ OTri P A → (s ∈ τ A ↔ hullOTriHom P φ hφ s ∈ τ B)

/-! ## ★(a) の後半 —— 分裂 `𝒪^×(A) × τ(A) ≅ 𝒪^▷(A)`

原文 (FrdI p.47):
> which is functorial in A; (b) for every isotropic hull A

★**原文は「hence determines」** —— 全単射性から**導かれる**と言っている。
★導出の要は `Proposition 2.2, (iii)` の単射性(`otri_div_eq_iff`)である。
-/

section Splitting

variable {τ : ∀ A : C, Submonoid (End A)} (hτ : IsCharacteristicSplitting P F τ)
  {A : C} (hA : IsIsotropic P A)

include hτ hA

/-- ★★★**Definition 2.3, (a) の分裂** —— `𝒪^×(A) × τ(A) → 𝒪^▷(A)`、
`(u, t) ↦ u · t` は**全単射**である。

★**全射性**は `charBij` の存在部分 ＋ `otri_div_eq_iff`(単射性)、
★**単射性**は `charBij` の一意性部分 ＋ `𝒞` の totally epimorphicity。 -/
theorem charSplitting_bijective :
    Function.Bijective
      (fun p : OTimes P A × τ A =>
        (⟨((p.1 : End A)) * ((p.2 : End A)),
          mul_mem (OTimes_le_OTri P A p.1.2) (hτ.le_otri A p.2.2)⟩ : OTri P A)) := by
  -- 単元は `Div = 0`、`degFr = 1`、`Base = 𝟙`
  have hunit : ∀ u : OTimes P A, P.Div (((u : End A)) : A ⟶ A) = 0 := by
    intro u
    haveI : IsIso (((u : End A)) : A ⟶ A) := (CategoryTheory.isUnit_iff_isIso _).mp u.2.2
    exact isIsometric_of_isIso P _
  -- `u * t = t ≫ u` なので `Div (u * t) = Div t`
  have hdiv : ∀ (u : OTimes P A) (t : OTri P A),
      P.Div ((((u : End A)) * ((t : End A)) : End A) : A ⟶ A)
        = P.Div (((t : End A)) : A ⟶ A) := by
    intro u t
    show P.Div ((((t : End A)) : A ⟶ A) ≫ (((u : End A)) : A ⟶ A)) = _
    rw [P.Div_comp, hunit u,
      show P.Base (((t : End A)) : A ⟶ A) = 𝟙 _ from by
        have h : P.Base (((t : End A)) : A ⟶ A) = P.Base (𝟙 A) := t.2.1
        rwa [P.Base_id] at h,
      MonoidOn.map_id,
      show P.degFr (((u : End A)) : A ⟶ A) = 1 from u.2.1.2]
    simp
  constructor
  · rintro ⟨u₁, t₁⟩ ⟨u₂, t₂⟩ h
    have h' : (((u₁ : End A)) * ((t₁ : End A)) : End A)
        = (((u₂ : End A)) * ((t₂ : End A)) : End A) := congrArg Subtype.val h
    -- `Div` を取ると `t₁` と `t₂` の `Div` が一致する
    have hd : P.Div (((t₁ : End A)) : A ⟶ A) = P.Div (((t₂ : End A)) : A ⟶ A) := by
      rw [← hdiv u₁ ⟨(t₁ : End A), hτ.le_otri A t₁.2⟩,
        ← hdiv u₂ ⟨(t₂ : End A), hτ.le_otri A t₂.2⟩]
      exact congrArg (fun z : End A => P.Div (z : A ⟶ A)) h'
    -- `charBij` の一意性で `t₁ = t₂`
    obtain ⟨w, -, huniq⟩ := hτ.charBij A hA ⟨(t₁ : End A), hτ.le_otri A t₁.2⟩
    have ht : t₁ = t₂ := (huniq t₁ rfl).trans (huniq t₂ hd.symm).symm
    subst ht
    -- 残りは epi による消去
    refine Prod.ext ?_ rfl
    refine Subtype.ext ?_
    haveI : Epi ((((t₁ : End A))) : A ⟶ A) := P.totEpiC _ _ _
    exact (cancel_epi ((((t₁ : End A))) : A ⟶ A)).mp h'
  · intro x
    obtain ⟨t, ht, -⟩ := hτ.charBij A hA x
    obtain ⟨u, hu⟩ := (otri_div_eq_iff P F hA x ⟨(t : End A), hτ.le_otri A t.2⟩).mp ht.symm
    exact ⟨⟨u, t⟩, Subtype.ext hu.symm⟩

omit hτ hA in
/-- ★**分裂の関手性 1** —— `otriLin` は `𝒪^×` を `𝒪^×` へ写す。

★`𝒪^▷(B)` の単元は `𝒪^▷(B)` の中で可逆(逆射も base-identity・linear)なので、
モノイド準同型 `otriLin` はそれを単元へ送る。 -/
theorem otriLin_otimes_mem {A B : C} (hA' : IsIsotropic P A) {φ : A ⟶ B}
    (hl : IsLinear P φ) (u : OTri P B) (hu : ((u : End B)) ∈ OTimes P B) :
    ((otriLin P F hA' hl u : End A)) ∈ OTimes P A := by
  haveI : IsIso ((((u : End B))) : B ⟶ B) := (CategoryTheory.isUnit_iff_isIso _).mp hu.2
  have hbu : P.Base ((((u : End B))) : B ⟶ B) = 𝟙 _ := by
    have h : P.Base ((((u : End B))) : B ⟶ B) = P.Base (𝟙 B) := u.2.1
    rwa [P.Base_id] at h
  have hinv : (inv ((((u : End B))) : B ⟶ B)) ∈ OTri P B := by
    refine ⟨?_, degFr_of_isIso P _⟩
    show P.Base (inv ((((u : End B))) : B ⟶ B)) = P.Base (𝟙 B)
    have h : P.Base (inv ((((u : End B))) : B ⟶ B)) ≫ P.Base ((((u : End B))) : B ⟶ B)
        = P.Base (𝟙 B) := by rw [← P.Base_comp, IsIso.inv_hom_id]
    rwa [hbu, Category.comp_id] at h
  -- `𝒪^▷(B)` の中での逆元
  have hmul : u * (⟨inv ((((u : End B))) : B ⟶ B), hinv⟩ : OTri P B) = 1 :=
    Subtype.ext (by show inv ((((u : End B))) : B ⟶ B) ≫ _ = _; simp)
  have hmul' : (⟨inv ((((u : End B))) : B ⟶ B), hinv⟩ : OTri P B) * u = 1 :=
    Subtype.ext (by show ((((u : End B))) : B ⟶ B) ≫ _ = _; simp)
  refine ⟨(otriLin P F hA' hl u).2, (CategoryTheory.isUnit_iff_isIso _).mpr ?_⟩
  refine ⟨(otriLin P F hA' hl ⟨inv ((((u : End B))) : B ⟶ B), hinv⟩ : End A), ?_, ?_⟩
  · have h := congrArg (fun z : OTri P A => ((z : End A)) )
      ((otriLin P F hA' hl).map_mul ⟨inv ((((u : End B))) : B ⟶ B), hinv⟩ u).symm
    simp only [hmul', map_one] at h
    exact h
  · have h := congrArg (fun z : OTri P A => ((z : End A)) )
      ((otriLin P F hA' hl).map_mul u ⟨inv ((((u : End B))) : B ⟶ B), hinv⟩).symm
    simp only [hmul, map_one] at h
    exact h

include hτ in
/-- ★★**分裂の関手性 2**(原文の「which is functorial in A」) ——
`otriLin` は `𝒪^×(B) × τ(B) → 𝒪^▷(B)` の分裂を
`𝒪^×(A) × τ(A) → 𝒪^▷(A)` の分裂へ写す。

★**`otriLin` がモノイド準同型であること**と、
★**`𝒪^×` と `τ` のどちらも `otriLin` で保たれること**の合わせ技である。 -/
theorem charSplitting_functorial {B : C} (hA' : IsIsotropic P A) {φ : A ⟶ B}
    (hl : IsLinear P φ) (u : OTri P B) (hu : ((u : End B)) ∈ OTimes P B)
    (t : OTri P B) (ht : ((t : End B)) ∈ τ B) :
    ((otriLin P F hA' hl (u * t) : End A))
        = ((otriLin P F hA' hl u : End A)) * ((otriLin P F hA' hl t : End A))
      ∧ ((otriLin P F hA' hl u : End A)) ∈ OTimes P A
      ∧ ((otriLin P F hA' hl t : End A)) ∈ τ A :=
  ⟨congrArg (fun z : OTri P A => ((z : End A))) ((otriLin P F hA' hl).map_mul u t),
   otriLin_otimes_mem P F hA' hl u hu,
   hτ.map_mem hA' hl t ht⟩

end Splitting

/-! ## ★非空虚性 —— unit-trivial かつすべてが isotropic なら `τ := 𝒪^▷` が取れる

★★**定義が空でないことを機械検証しておく。**
★`𝒪^×(A) = ⊥` なら `𝒪^▷(A)^char = 𝒪^▷(A)` なので、`τ := 𝒪^▷` が (a) を満たす。
★すべての対象が isotropic なら isotropic hull は同型なので (b) も自明に満たす。
-/

/-- ★**isotropic な対象の isotropic hull は同型**。

`IsIsotropicHull` は isometric な pre-step だから、`IsIsotropic` の定義そのもの。 -/
theorem isIso_of_isotropicHull_of_isotropic {A B : C} {φ : A ⟶ B}
    (hA : IsIsotropic P A) (hφ : IsIsotropicHull P φ) : IsIso φ :=
  hA B φ hφ.1 hφ.2.1

include F in
/-- ★★★**Definition 2.3 の非空虚性** —— すべての対象が isotropic かつ unit-trivial な
Frobenioid では、`τ := 𝒪^▷(−)` 自身が characteristic splitting である。 -/
theorem isCharacteristicSplitting_otri
    (hiso : ∀ A : C, IsIsotropic P A) (hut : ∀ A : C, IsUnitTrivial P A) :
    IsCharacteristicSplitting P F (OTri P) where
  le_otri _ := le_rfl
  map_mem := by
    intro A B hA φ hl t _
    exact (otriLin P F hA hl t).2
  charBij A _ x := by
    refine ⟨x, rfl, ?_⟩
    rintro ⟨t, htm⟩ (ht : P.Div _ = P.Div _)
    obtain ⟨u, hu⟩ := (otri_div_eq_iff P F (hiso A) ⟨t, htm⟩ x).mp ht
    have hu1 : ((u : End A)) = 1 := by
      have hb : (u : End A) ∈ (⊥ : Submonoid (End A)) := by rw [← hut A]; exact u.2
      simpa using hb
    refine Subtype.ext ?_
    have hu' : t = ((x : End A)) ≫ ((u : End A)) := hu
    show t = (x : End A)
    rw [hu', hu1]
    simp
  hullMem {A B φ} hφ t htm := by
    haveI : IsIso φ := isIso_of_isotropicHull_of_isotropic P (hiso A) hφ
    -- `s := φ ≫ t ≫ φ⁻¹`
    have hsq : (φ ≫ t ≫ inv φ) ≫ φ = φ ≫ t := by simp
    have hbt : P.Base t = 𝟙 _ := by
      have h : P.Base t = P.Base (𝟙 B) := htm.1
      rwa [P.Base_id] at h
    have hdt : P.degFr t = 1 := htm.2
    have hbs : IsBaseIdentity P (φ ≫ t ≫ inv φ) := by
      show P.Base (φ ≫ t ≫ inv φ) = P.Base (𝟙 A)
      rw [P.Base_comp, P.Base_comp, hbt, Category.id_comp, ← P.Base_comp,
        IsIso.hom_inv_id]
    have hds : P.degFr (φ ≫ t ≫ inv φ) = 1 := by
      rw [P.degFr_comp, P.degFr_comp, hdt, degFr_of_isIso P φ,
        degFr_of_isIso P (inv φ)]
      simp
    exact ⟨⟨φ ≫ t ≫ inv φ, hbs, hds⟩,
      (hullOTriMap_uniq P φ hφ (φ ≫ t ≫ inv φ) t hsq).symm⟩
  hullPullback {A B φ} hφ s hs :=
    ⟨fun _ => hullOTriHom_mem P φ hφ s hs, fun _ => hs⟩

/-! ## ★★★出典の紐付け(`.src`) -/

/-- ★★★**[FrdI] Definition 2.3** —— 4 主張すべてが実装された。

| # | 主張 | 実装 |
|---|---|---|
| 1 | 部分関手 | `IsCharacteristicSplitting.le_otri` / `.map_mem` |
| 2 | (a) `τ(A) ≅ 𝒪^▷(A)^char` | `IsCharacteristicSplitting.charBij` |
| 3 | (a) 分裂 `𝒪^×(A) × τ(A) ≅ 𝒪^▷(A)` | `charSplitting_bijective` |
| 4 | (b) isotropic hull の像 | `IsCharacteristicSplitting.hullMem` |

★非空虚性は `isCharacteristicSplitting_otri` で機械検証した。 -/
def IsCharacteristicSplitting.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 47, item := "Definition 2.3",
    sectionId := "frdi-def-2-3" }

end ABC3.Found.FrdI
