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

/-! ## ★7. 辞書の中心 —— **co-angular pre-step は `𝒞^birat` で同型になる**

原文 (FrdI p.83):
> morphism of a given Frobenius degree; isometry; pre-step; base-isomorphism) of

★★**逆射は「添字を `φ` 自身に取って恒等射を置く」だけ**である ——
`Hom^birat(B, A)` の添字は `B` へ入る co-angular pre-step なので、
★`φ : A ⟶ B` 自身が添字になり、そこでの射 `𝟙_A` が逆射を与える。

★★**2 つの合成則の計算はどちらも同じ形**で落ちる:
1. `biratPull_sq` の四角形と **`φ`(または `𝟙`)が mono** であることから `γ = α`
2. `HomBirat.mk_map`(前合成で移しても同じ元)で `𝟙` の代表に戻す -/

theorem birat_inv_comp {A B : C} (φ : A ⟶ B) (hc : IsCoAngular P φ) (hs : IsPreStep P φ) :
    compBirat P G G.core (toHomBirat (P := P) (G := G) φ)
        (HomBirat.mk (idxBiratMk P G φ hc hs) (𝟙 A))
      = toHomBirat (P := P) (G := G) (𝟙 A) := by
  haveI : Mono φ := G.core.preStepMono φ hs
  have hsq := biratPull_sq G.core (idxBiratOne P G A) φ (idxBiratMk P G φ hc hs)
  have hga : biratPullGamma G.core (idxBiratOne P G A) φ (idxBiratMk P G φ hc hs)
      = biratPullAlpha G.core (idxBiratOne P G A) φ (idxBiratMk P G φ hc hs) :=
    (cancel_mono φ).mp hsq
  rw [show toHomBirat (P := P) (G := G) φ = HomBirat.mk (idxBiratOne P G A) φ from rfl,
    compBirat_mk G.core (idxBiratOne P G A) φ (idxBiratMk P G φ hc hs)]
  refine Eq.trans ?_ (HomBirat.mk_map
    (idxBiratHomMk (Z := idxBiratOne P G A)
      (W := biratPullIdx G.core (idxBiratOne P G A) φ (idxBiratMk P G φ hc hs))
      (biratPullGamma G.core (idxBiratOne P G A) φ (idxBiratMk P G φ hc hs))
      (biratPullGamma_coAngular G.core (idxBiratOne P G A) φ (idxBiratMk P G φ hc hs))
      (biratPullGamma_preStep G.core (idxBiratOne P G A) φ (idxBiratMk P G φ hc hs))
      rfl) (𝟙 A))
  refine congrArg (HomBirat.mk
    (biratPullIdx G.core (idxBiratOne P G A) φ (idxBiratMk P G φ hc hs))) ?_
  rw [idxBiratHomMk_left]
  exact congrArg (fun t => t ≫ 𝟙 A) hga.symm

theorem birat_comp_inv {A B : C} (φ : A ⟶ B) (hc : IsCoAngular P φ) (hs : IsPreStep P φ) :
    compBirat P G G.core (HomBirat.mk (idxBiratMk P G φ hc hs) (𝟙 A))
        (toHomBirat (P := P) (G := G) φ)
      = toHomBirat (P := P) (G := G) (𝟙 B) := by
  have hsq := biratPull_sq G.core (idxBiratMk P G φ hc hs) (𝟙 A) (idxBiratOne P G A)
  have hga : biratPullGamma G.core (idxBiratMk P G φ hc hs) (𝟙 A) (idxBiratOne P G A)
      = biratPullAlpha G.core (idxBiratMk P G φ hc hs) (𝟙 A) (idxBiratOne P G A) :=
    (Category.comp_id _).symm.trans (hsq.trans (Category.comp_id _))
  have hcc : IsCoAngular P
      (biratPullGamma G.core (idxBiratMk P G φ hc hs) (𝟙 A) (idxBiratOne P G A) ≫ φ) :=
    G.core.coAngularComp _ _
      (biratPullGamma_coAngular G.core (idxBiratMk P G φ hc hs) (𝟙 A) (idxBiratOne P G A)) hc
  have hcs : IsPreStep P
      (biratPullGamma G.core (idxBiratMk P G φ hc hs) (𝟙 A) (idxBiratOne P G A) ≫ φ) :=
    IsPreStep.comp P
      (biratPullGamma_preStep G.core (idxBiratMk P G φ hc hs) (𝟙 A) (idxBiratOne P G A)) hs
  rw [show toHomBirat (P := P) (G := G) φ = HomBirat.mk (idxBiratOne P G A) φ from rfl,
    compBirat_mk G.core (idxBiratMk P G φ hc hs) (𝟙 A) (idxBiratOne P G A)]
  refine Eq.trans ?_ (HomBirat.mk_map
    (idxBiratHomMk (Z := idxBiratOne P G B)
      (W := biratPullIdx G.core (idxBiratMk P G φ hc hs) (𝟙 A) (idxBiratOne P G A))
      (biratPullGamma G.core (idxBiratMk P G φ hc hs) (𝟙 A) (idxBiratOne P G A) ≫ φ)
      hcc hcs (Category.comp_id _)) (𝟙 B))
  refine congrArg (HomBirat.mk
    (biratPullIdx G.core (idxBiratMk P G φ hc hs) (𝟙 A) (idxBiratOne P G A))) ?_
  rw [idxBiratHomMk_left]
  exact (congrArg (fun t => t ≫ φ) hga.symm).trans (Category.comp_id _).symm

/-- ★★★**[FrdI] Proposition 4.4, (iv) の中心** —— `𝒞` の co-angular pre-step は
`𝒞^birat` の**同型**になる。

原文 (FrdI p.83):
> morphism of a given Frobenius degree; isometry; pre-step; base-isomorphism) of

★★**これが「`𝒞^birat` は co-angular pre-step を反転して作った圏である」**
という構成の意味そのものである。 -/
theorem birat_isIso_of_coaPre {A B : C} (φ : A ⟶ B) (hc : IsCoAngular P φ)
    (hs : IsPreStep P φ) : IsIso ((toBiratCat P G).map φ) :=
  ⟨HomBirat.mk (idxBiratMk P G φ hc hs) (𝟙 A),
    birat_inv_comp φ hc hs, birat_comp_inv φ hc hs⟩

/-! ## ★8. 残り —— 辞書の 4 項と Frobenioid 構造

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


/-! ## ★9. ★★★★[FrdI] Proposition 4.8, (i) —— isotropic 型は `𝒞^birat` へ渡る

原文 (FrdI p.88):
> Assertion (i) follows formally from Proposition 4.4, (iv). To prove assertion

★★**`𝒞^birat` ではすべての射が等長**(零因子が 1 元単系)なので、
**isotropic 性は「pre-step がすべて同型」**に潰れる。
★`𝒞^birat` の pre-step は代表元 `(a, φ)` で `φ` が `𝒞` の pre-step であること。
★`𝒞` が isotropic 型なら `Proposition 1.4, (i)` で `φ` は co-angular、
したがって `𝒞^birat` で同型(`birat_isIso_of_coaPre`)。
★★あとは `a` も同型なので、`f = [a]⁻¹ ≫ [φ]` が同型になる。 -/

/-- ★★**代表元は `a` と `φ` の合成に分かれる** —— `[a] ≫ [a⁻¹ ≫ φ] = [φ]`。 -/
theorem birat_toHom_comp_mk {A B A' : C} (a : A' ⟶ A) (hac : IsCoAngular P a)
    (has : IsPreStep P a) (φ : A' ⟶ B) :
    compBirat P G G.core (toHomBirat (P := P) (G := G) a)
        (HomBirat.mk (idxBiratMk P G a hac has) φ)
      = toHomBirat (P := P) (G := G) φ := by
  haveI : Mono a := G.core.preStepMono a has
  have hsq := biratPull_sq G.core (idxBiratOne P G A') a (idxBiratMk P G a hac has)
  have hga : biratPullGamma G.core (idxBiratOne P G A') a (idxBiratMk P G a hac has)
      = biratPullAlpha G.core (idxBiratOne P G A') a (idxBiratMk P G a hac has) :=
    (cancel_mono a).mp hsq
  rw [show toHomBirat (P := P) (G := G) a = HomBirat.mk (idxBiratOne P G A') a from rfl,
    compBirat_mk G.core (idxBiratOne P G A') a (idxBiratMk P G a hac has)]
  refine Eq.trans ?_ (HomBirat.mk_map
    (idxBiratHomMk (Z := idxBiratOne P G A')
      (W := biratPullIdx G.core (idxBiratOne P G A') a (idxBiratMk P G a hac has))
      (biratPullGamma G.core (idxBiratOne P G A') a (idxBiratMk P G a hac has))
      (biratPullGamma_coAngular G.core (idxBiratOne P G A') a (idxBiratMk P G a hac has))
      (biratPullGamma_preStep G.core (idxBiratOne P G A') a (idxBiratMk P G a hac has))
      rfl) φ)
  refine congrArg (HomBirat.mk
    (biratPullIdx G.core (idxBiratOne P G A') a (idxBiratMk P G a hac has))) ?_
  rw [idxBiratHomMk_left]
  exact congrArg (fun t => t ≫ φ) hga.symm

/-- ★★**`𝒞^birat` の pre-step は代表元の `φ` が `𝒞` の pre-step であること**。 -/
theorem birat_preStep_rep {A B : C} (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B)
    (h : IsPreStep (biratPre P G) (HomBirat.mk Z φ)) : IsPreStep P φ := by
  haveI hZ : IsIso (P.Base Z.unop.hom.hom) := Z.unop.hom.property.2.2
  refine ⟨?_, ?_⟩
  · have hd : biratDeg (HomBirat.mk Z φ) = 1 := h.1
    rw [biratDeg_mk] at hd
    exact hd
  · have hb : IsIso (biratBase (HomBirat.mk Z φ)) := h.2
    rw [biratBase_mk, sliceBaseOf_eq] at hb
    haveI := hb
    have : IsIso (P.Base Z.unop.hom.hom ≫ (inv (P.Base Z.unop.hom.hom) ≫ P.Base φ)) :=
      IsIso.comp_isIso' inferInstance hb
    rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp] at this
    exact this

include P in
/-- ★★★★**[FrdI] Proposition 4.8, (i)** —— `𝒞` が isotropic 型なら `𝒞^birat` もそう。

原文 (FrdI p.88):
> Assertion (i) follows formally from Proposition 4.4, (iv). To prove assertion -/
theorem birat_isOfIsotropicType (hiso : ∀ X : C, IsIsotropic P X)
    (X : BiratCat P G) : IsIsotropic (biratPre P G) X := by
  intro Y f _ hs
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep f
  -- ★代表元の `φ` は `𝒞` の pre-step、しかも co-angular
  have hφs : IsPreStep P φ := birat_preStep_rep Z φ hs
  have hφc : IsCoAngular P φ := prop_1_4_i P φ (fun W _ => hiso W)
  haveI : IsIso ((toBiratCat P G).map φ) := birat_isIso_of_coaPre φ hφc hφs
  -- ★添字の構造射 `a` も co-angular pre-step なので同型
  haveI : IsIso ((toBiratCat P G).map Z.unop.hom.hom) :=
    birat_isIso_of_coaPre Z.unop.hom.hom Z.unop.hom.property.1 Z.unop.hom.property.2
  -- ★`[a] ≫ f = [φ]` なので `f` は同型
  have hcomp : (toBiratCat P G).map Z.unop.hom.hom ≫ HomBirat.mk Z φ
      = (toBiratCat P G).map φ := by
    refine Eq.trans ?_ (birat_toHom_comp_mk Z.unop.hom.hom Z.unop.hom.property.1
      Z.unop.hom.property.2 φ)
    rfl
  haveI : IsIso ((toBiratCat P G).map Z.unop.hom.hom ≫ HomBirat.mk Z φ) := by
    rw [hcomp]
    exact birat_isIso_of_coaPre φ hφc hφs
  exact IsIso.of_isIso_comp_left ((toBiratCat P G).map Z.unop.hom.hom) (HomBirat.mk Z φ)

/-! ## ★10. ★★★★[FrdI] Proposition 4.4, (ii) —— `𝒞 → 𝒞^birat` は忠実

原文 (FrdI p.84):
> C →Cbirat is faithful and determines an injection O▷(A)gp →O×(Abirat), for

★★**帰納極限の等号判定 ＋ `𝒞` が totally epimorphic**、それだけである。 -/

include P in
theorem toBiratCat_faithful : (toBiratCat P G).Faithful where
  map_injective {A B} {f g} h := by
    obtain ⟨V, u, hu⟩ := HomBirat.eq_iff_same.mp h
    haveI : Epi u.unop.left.hom := P.totEpiC _ _ _
    exact (cancel_epi u.unop.left.hom).mp hu

/-! ## ★11. ★★★★辞書の「同型 ⟹ co-angular pre-step」

原文 (FrdI p.83):
> morphism of a given Frobenius degree; isometry; pre-step; base-isomorphism) of

★★**筋**(2026-08-17、紙の上で確定):
1. 逆射の代表を `(b, χ)` と書くと、`[b] ≫ (逆射) = [χ]` なので `[χ ≫ φ] = [b]`。
   ★**`𝒞 → 𝒞^birat` が忠実**(上)なので `χ ≫ φ = b`。
2. `φ` を `Definition 1.3, (v), (b)` で `φ = β ≫ α`(`β` co-angular pre-step、
   `α` **等長** pre-step)と分解する。
3. ★★★**`b` の co-angular 性を分解 `b = (χ ≫ β) ≫ α ≫ 𝟙` に当てる** ——
   真ん中が等長 pre-step、最後が線型かつ底同型なので、`α` は同型。
4. よって `φ = β ≫ (同型)` は co-angular pre-step。

★★**要点は 3** —— 「co-angular の定義は**真ん中の等長 pre-step が同型**」なので、
`b` の co-angular 性が `φ` の等長部分をそのまま潰す。 -/

include P in
/-- ★★★**同型なら代表元は co-angular pre-step** —— 辞書の「同型」の `⟹`。 -/
theorem birat_coaPre_of_isIso {A B : C} (φ : A ⟶ B)
    (hiso : IsIso ((toBiratCat P G).map φ)) :
    IsCoAngular P φ ∧ IsPreStep P φ := by
  haveI := hiso
  haveI := toBiratCat_faithful (P := P) (G := G)
  -- ★段 0: `φ` は pre-step
  have hdeg : P.degFr φ = 1 := by
    have h := biratDeg_toHomBirat (P := P) (G := G) φ
    have h1 : biratDeg ((toBiratCat P G).map φ) = 1 := by
      have := degFr_of_isIso (biratPre P G) ((toBiratCat P G).map φ)
      exact this
    rw [← h]; exact h1
  have hbase : IsBaseIsomorphism P φ := by
    have h1 : IsIso ((biratPre P G).Base ((toBiratCat P G).map φ)) :=
      isBaseIsomorphism_of_isIso (biratPre P G) ((toBiratCat P G).map φ)
    show IsIso (P.Base φ)
    rw [← biratBase_toHomBirat (P := P) (G := G) φ]
    exact h1
  have hφs : IsPreStep P φ := ⟨hdeg, hbase⟩
  -- ★段 1: 逆射の代表から `χ ≫ φ = b`
  obtain ⟨g, hg1, hg2⟩ := hiso.out
  obtain ⟨W, χ, hWχ⟩ := HomBirat.exists_rep (P := P) (G := G) g
  have hb : χ ≫ φ = W.unop.hom.hom := by
    refine (toBiratCat P G).map_injective ?_
    have hstep : (toBiratCat P G).map W.unop.hom.hom ≫ g = (toBiratCat P G).map χ := by
      refine Eq.trans ?_ (birat_toHom_comp_mk W.unop.hom.hom W.unop.hom.property.1
        W.unop.hom.property.2 χ)
      rw [← hWχ]
      rfl
    have h2 : (toBiratCat P G).map χ ≫ (toBiratCat P G).map φ
        = (toBiratCat P G).map W.unop.hom.hom :=
      calc (toBiratCat P G).map χ ≫ (toBiratCat P G).map φ
          = ((toBiratCat P G).map W.unop.hom.hom ≫ g) ≫ (toBiratCat P G).map φ :=
            congrArg (fun t => t ≫ (toBiratCat P G).map φ) hstep.symm
        _ = (toBiratCat P G).map W.unop.hom.hom ≫ (g ≫ (toBiratCat P G).map φ) :=
            Category.assoc _ _ _
        _ = (toBiratCat P G).map W.unop.hom.hom ≫ 𝟙 _ :=
            congrArg (fun t => (toBiratCat P G).map W.unop.hom.hom ≫ t) hg2
        _ = (toBiratCat P G).map W.unop.hom.hom := Category.comp_id _
    rw [← (toBiratCat P G).map_comp] at h2
    exact h2
  -- ★段 2: `φ` を分解する
  obtain ⟨X, β, α, heq, hβc, hβs, hαi, hαs⟩ := G.core.preStepFactor φ hφs
  -- ★段 3: `b` の co-angular 性で `α` は同型
  have hbc : IsCoAngular P W.unop.hom.hom := W.unop.hom.property.1
  haveI hαiso : IsIso α := by
    refine hbc X B (χ ≫ β) α (𝟙 B) ?_ ?_ hαi hαs ?_
    · calc W.unop.hom.hom = χ ≫ φ := hb.symm
        _ = χ ≫ (β ≫ α) := congrArg (fun t => χ ≫ t) heq
        _ = (χ ≫ β) ≫ α := (Category.assoc _ _ _).symm
        _ = (χ ≫ β) ≫ α ≫ 𝟙 B := by rw [Category.comp_id]
    · exact isLinear_of_isIso P (𝟙 B)
    · exact Or.inl (isBaseIsomorphism_of_isIso P (𝟙 B))
  -- ★段 4: `φ = β ≫ (同型)`
  refine ⟨?_, hφs⟩
  rw [heq]
  exact IsCoAngular.comp P G.core hβc (isCoAngular_of_isIso P α)

include P in
/-- ★★★★**[FrdI] Proposition 4.4, (iv) の「同型」の条** ——
`𝒞` の射が `𝒞^birat` で同型になるのは、ちょうど **co-angular pre-step** のとき。

原文 (FrdI p.83):
> morphism of a given Frobenius degree; isometry; pre-step; base-isomorphism) of -/
theorem birat_isIso_iff {A B : C} (φ : A ⟶ B) :
    IsIso ((toBiratCat P G).map φ) ↔ (IsCoAngular P φ ∧ IsPreStep P φ) :=
  ⟨birat_coaPre_of_isIso φ, fun h => birat_isIso_of_coaPre φ h.1 h.2⟩

end ABC3.Found.FrdI
