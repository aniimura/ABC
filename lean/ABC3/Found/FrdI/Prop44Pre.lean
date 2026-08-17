import ABC3.Found.FrdI.Prop44

/-!
# [FrdI] Proposition 4.4, (ii) の土台 —— `𝒞^birat` の pre-Frobenioid 構造

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.83。

原文 (FrdI p.83):
> (ii) The functor Cbirat →F0D of (i) determines a structure of Frobenioid

★`Prop44.lean` で圏 `𝒞^birat` と関手 `𝒞^birat → 𝔽_{0_𝒟}` は取れている。
★★**ここでは `PreFrobenioid` の 6 フィールドを埋める。**

## ★★6 フィールドの内訳(測定、2026-08-17)

| フィールド | 手 |
|---|---|
| `toElem` | `biratToElem`(取得済み) |
| `divisorial` | ★**台が 1 元**なので 4 条件すべて自明 |
| `totEpiD` / `connectedD` | `𝒟` は変わらないので `P` のものがそのまま |
| `connectedC` | `toBiratCat` は**対象について恒等**なので zigzag を送るだけ |
| `totEpiC` | ★★**唯一の中身**(下) |

## ★★★`𝒞^birat` が totally epimorphic であること

★★**筋**:
1. `g`・`h` を**同じ添字**の代表元に揃える(添字圏が filtered)
2. 合成 `f ≫ g`・`f ≫ h` は**同じ引き戻し**(`biratPullIdx` / `biratPullAlpha`)を使うので、
   等式は `mk V (α ≫ ψ) = mk V (α ≫ ψ′)` になる
3. 帰納極限の等号判定(`eq_iff_same`)で**ある遷移射 `u` を当てた後で一致**
4. ★**`𝒞` が totally epimorphic なので `u ≫ α` は epi**、よって `ψ = ψ′`

★★**要点は 2**である —— `f` の代表を固定すれば、`g` と `h` の合成が
**同じ添字・同じ `α`** で書けるので、話が `𝒞` の中の等式に落ちる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory CategoryTheory.Limits

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (G : Frobenioid P)

/-! ## ★1. `0_𝒟` は divisorial —— 台が 1 元だから -/

/-- ★**1 元単系は divisorial**。★4 条件(integral・saturated・characteristic 型・sharp)が
すべて `Subsingleton.elim` で落ちる。 -/
theorem subsingleton_gp_punit : Subsingleton (Gp (PUnit.{w + 1})) := by
  refine ⟨fun x y => ?_⟩
  induction x using AddLocalization.induction_on with
  | H a =>
    induction y using AddLocalization.induction_on with
    | H b =>
      exact congrArg₂ AddLocalization.mk (Subsingleton.elim a.1 b.1)
        (Subsingleton.elim a.2 b.2)

theorem isDivisorial_punit : IsDivisorial (PUnit.{w + 1}) := by
  haveI := subsingleton_gp_punit.{w}
  refine ⟨⟨?_, ?_, ?_⟩, ?_⟩
  · intro a b _
    exact Subsingleton.elim a b
  · intro a _ _ _
    exact ⟨(PUnit.unit : PUnit.{w + 1}), Subsingleton.elim _ a⟩
  · intro a b _
    exact ⟨0, ⟨isAddUnit_zero, Subsingleton.elim _ _⟩,
      fun _ _ => Subsingleton.elim (α := PUnit.{w + 1}) _ _⟩
  · intro a _
    exact Subsingleton.elim a 0

variable (D) in
theorem trivialOn_isDivisorialOn : (trivialOn.{v, u, w} D).IsDivisorialOn :=
  fun _ => isDivisorial_punit

/-! ## ★2. 代表元を揃える -/

variable {P G}

/-- ★**2 つの射を同じ添字の代表元に揃える** —— 添字圏が filtered だから。 -/
theorem HomBirat.exists_rep_pair {B E : C} (g h : HomBirat P G B E) :
    ∃ (W : IdxBirat P G B) (ψ ψ' : W.unop.left.obj ⟶ E),
      HomBirat.mk W ψ = g ∧ HomBirat.mk W ψ' = h := by
  obtain ⟨W₁, ψ₁, h₁⟩ := g.exists_rep
  obtain ⟨W₂, ψ₂, h₂⟩ := h.exists_rep
  refine ⟨IsFiltered.max W₁ W₂,
    (IsFiltered.leftToMax W₁ W₂).unop.left.hom ≫ ψ₁,
    (IsFiltered.rightToMax W₁ W₂).unop.left.hom ≫ ψ₂, ?_, ?_⟩
  · exact (HomBirat.mk_map (IsFiltered.leftToMax W₁ W₂) ψ₁).trans h₁
  · exact (HomBirat.mk_map (IsFiltered.rightToMax W₁ W₂) ψ₂).trans h₂

/-- ★**同じ添字での等号判定**。 -/
theorem HomBirat.eq_iff_same {A B : C} {Z : IdxBirat P G A} {φ ψ : Z.unop.left.obj ⟶ B} :
    HomBirat.mk Z φ = HomBirat.mk Z ψ ↔
      ∃ (V : IdxBirat P G A) (u : Z ⟶ V), u.unop.left.hom ≫ φ = u.unop.left.hom ≫ ψ := by
  refine (HomColim.eq_iff_same (homFunctorBirat P G A B)).trans ⟨?_, ?_⟩
  · rintro ⟨V, u, h⟩
    exact ⟨V, u, ULift.up_injective h⟩
  · rintro ⟨V, u, h⟩
    exact ⟨V, u, congrArg ULift.up h⟩

/-! ## ★3. `𝒞^birat` は totally epimorphic -/

variable (P G)

/-- ★★★**`𝒞^birat` は totally epimorphic**。

★合成の定義(`compBirat_mk`)が **`f` の代表を固定すれば同じ引き戻しを使う**ので、
等式が `𝒞` の中の等式に落ち、`𝒞` の totally epimorphic 性で閉じる。 -/
theorem birat_totEpi : IsTotallyEpimorphic (BiratCat P G) := by
  intro A B f
  refine ⟨fun {E} g h hgh => ?_⟩
  obtain ⟨Z, φ, hφ⟩ := HomBirat.exists_rep (P := P) (G := G) f
  obtain ⟨W, ψ, ψ', hψ, hψ'⟩ := HomBirat.exists_rep_pair
    (P := P) (G := G) (B := biratDown P G B) (E := biratDown P G E) g h
  -- ★合成を代表元で書く
  have hc : ∀ χ : W.unop.left.obj ⟶ biratDown P G E,
      compBirat P G G.core (HomBirat.mk Z φ) (HomBirat.mk W χ)
        = HomBirat.mk (biratPullIdx G.core Z φ W) (biratPullAlpha G.core Z φ W ≫ χ) :=
    fun χ => compBirat_mk G.core Z φ W χ
  have heq : HomBirat.mk (biratPullIdx G.core Z φ W) (biratPullAlpha G.core Z φ W ≫ ψ)
      = HomBirat.mk (biratPullIdx G.core Z φ W) (biratPullAlpha G.core Z φ W ≫ ψ') := by
    rw [← hc ψ, ← hc ψ', hφ, hψ, hψ']
    exact hgh
  obtain ⟨V, u, hu⟩ := HomBirat.eq_iff_same.mp heq
  -- ★`u ≫ α` は epi
  haveI : Epi (u.unop.left.hom ≫ biratPullAlpha G.core Z φ W) := P.totEpiC _ _ _
  refine hψ ▸ hψ' ▸ congrArg (HomBirat.mk W) ?_
  refine (cancel_epi (u.unop.left.hom ≫ biratPullAlpha G.core Z φ W)).mp ?_
  rw [Category.assoc, Category.assoc]
  exact hu

/-! ## ★4. `𝒞^birat` は connected —— 対象は `𝒞` と同じ -/

theorem birat_isConnected : IsConnected (BiratCat P G) := by
  haveI := P.connectedC
  haveI : Nonempty (BiratCat P G) := (inferInstance : Nonempty C)
  refine zigzag_isConnected ?_
  intro X Y
  exact zigzag_obj_of_zigzag (toBiratCat P G)
    (@isPreconnected_zigzag C _ _ (biratDown P G X) (biratDown P G Y))

/-! ## ★5. ★★★★`𝒞^birat` の pre-Frobenioid 構造 -/

/-- ★★★★**[FrdI] Proposition 4.4, (ii) の第 1 段** ——
`𝒞^birat → 𝔽_{0_𝒟}` は `𝒞^birat` に **pre-Frobenioid の構造**を定める。

原文 (FrdI p.83):
> (ii) The functor Cbirat →F0D of (i) determines a structure of Frobenioid

★6 フィールドのうち **4 つは `P` のものと 1 元単系の自明性で埋まり**、
実質は `totEpiC`(`birat_totEpi`)だけである。 -/
noncomputable def biratPre : PreFrobenioid (BiratCat P G) (trivialOn.{v, u, w} D) where
  toElem := biratToElem P G
  divisorial := trivialOn_isDivisorialOn D
  totEpiC := birat_totEpi P G
  totEpiD := P.totEpiD
  connectedC := birat_isConnected P G
  connectedD := P.connectedD

/-! ### ★不変量の計算則 -/

variable {P G}

@[simp] theorem biratPre_degFr {A B : BiratCat P G} (f : A ⟶ B) :
    (biratPre P G).degFr f = biratDeg f := rfl

@[simp] theorem biratPre_Base {A B : BiratCat P G} (f : A ⟶ B) :
    (biratPre P G).Base f = biratBase f := rfl

/-- ★★**`𝒞^birat` では零因子は情報を持たない** —— 台が 1 元だから。 -/
theorem biratPre_Div_eq_zero {A B : BiratCat P G} (f : A ⟶ B) :
    (biratPre P G).Div f = 0 :=
  Subsingleton.elim (α := PUnit.{w + 1}) _ _

/-- ★★★**`𝒞^birat` ではすべての射が等長** —— 上の帰結。

原文 (FrdI p.83):
> of a given Frobenius degree; isometry; pre-step; base-isomorphism) of

★★これが原文の辞書の「isometry ⟺ **arbitrary morphism**」の中身である。 -/
theorem birat_isIsometric {A B : BiratCat P G} (f : A ⟶ B) :
    IsIsometric (biratPre P G) f :=
  biratPre_Div_eq_zero f

/-! ## ★6. `Proposition 4.4, (iv)` の辞書のうち**不変量で決まる 3 項**

原文 (FrdI p.83):
> of a given Frobenius degree; isometry; pre-step; base-isomorphism) of

★★`𝒞 → 𝒞^birat` は**次数と底をそのまま保つ**(`biratDeg_toHomBirat` /
`biratBase_toHomBirat`)ので、この 3 項は**両辺が同じ式**になる。 -/

/-- ★**与えられた Frobenius 次数** —— そのまま渡る。 -/
theorem birat_degFr_iff {A B : C} (φ : A ⟶ B) (n : ℕ+) :
    (biratPre P G).degFr ((toBiratCat P G).map φ) = n ↔ P.degFr φ = n := by
  show biratDeg (toHomBirat (P := P) (G := G) φ) = n ↔ _
  rw [biratDeg_toHomBirat]

/-- ★**linear** —— そのまま渡る。 -/
theorem birat_isLinear_iff {A B : C} (φ : A ⟶ B) :
    IsLinear (biratPre P G) ((toBiratCat P G).map φ) ↔ IsLinear P φ :=
  birat_degFr_iff φ 1

/-- ★**base-isomorphism** —— そのまま渡る。 -/
theorem birat_isBaseIsomorphism_iff {A B : C} (φ : A ⟶ B) :
    IsBaseIsomorphism (biratPre P G) ((toBiratCat P G).map φ) ↔ IsBaseIsomorphism P φ := by
  show IsIso (biratBase (toHomBirat (P := P) (G := G) φ)) ↔ _
  rw [biratBase_toHomBirat]
  rfl

/-- ★★**pre-step** —— そのまま渡る(linear ∧ base-isomorphism)。

原文 (FrdI p.83):
> of a given Frobenius degree; isometry; pre-step; base-isomorphism) of
-/
theorem birat_isPreStep_iff {A B : C} (φ : A ⟶ B) :
    IsPreStep (biratPre P G) ((toBiratCat P G).map φ) ↔ IsPreStep P φ :=
  and_congr (birat_isLinear_iff φ) (birat_isBaseIsomorphism_iff φ)

/-! ## ★7. 残り —— 辞書の 5 項と Frobenioid 構造

★★**残っているのは、単に不変量で決まらない 5 項**である:

| `𝒞^birat` で | `𝒞` で |
|---|---|
| 同型 | ★**co-angular pre-step** |
| Frobenius 型 | ★**co-angular base-isomorphism** |
| pull-back | ★**co-angular linear** |
| co-angular | co-angular |
| base-identity 自己射 | ★(α, φ) の対で base-equivalent |

★★どれも **`𝒞^birat` の射が「co-angular pre-step で割った射」**であることを使う。
★その先に `Frobenioid` の 21 + 2 フィールド(原文の「routine exercise」)がある。
-/

end ABC3.Found.FrdI
