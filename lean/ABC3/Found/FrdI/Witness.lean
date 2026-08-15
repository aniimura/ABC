import ABC3.Found.FrdI.MorphismTypes

/-!
# ★[FrdI] Definition 1.1 / 1.2 の**非退化 witness**

定義を書いただけでは「型は付くが中身が空」を排除できない。ここでは

1. **本物の pre-Frobenioid を1つ組み立て**(`Definition 1.1, (iv)` の仮定を全部満たす)、
2. その上で `Definition 1.2` の各語について**満たす射と満たさない射の両方**を示す。

## 組み立てるもの

* 底の圏 `𝒟 := Vee` —— V 字の半順序(`left ≤ top`、`right ≤ top`、`left`/`right` は比較不能)。
  ★`CategoryVocabulary.lean` で **totally epimorphic** であることを示してある。
  ★**同型でない射を持つ**ので `base-isomorphism` の反例が取れる。
* divisor monoid `Φ :=` `ℕ` への定数関手。★`ℕ` が **divisorial** であることは
  `MonoidVocabulary.lean` の `isDivisorial_nat`(saturated を含む)で示してある。
* `𝒞 := 𝔽_Φ` 自身、`𝒞 → 𝔽_Φ` は恒等関手。

★この witness が退化していないことは、下の `linear` / `isometric` /
`base-isomorphism` / `pre-step` / `step` の**それぞれについて両方向**の例が
取れることで確かめる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

universe v u w

/-! ### 定数関手から `monoid on 𝒟` を作る -/

/-- 定数関手 `Φ(A) = M`、`α* = id` から作る `monoid on 𝒟`。

(a) `id` は characteristically injective、(b) `id` は全単射、なので2条件は自動。 -/
@[reducible]
def MonoidOn.const (D : Type u) [Category.{v} D] (M : Type w) [AddCommMonoid M] :
    MonoidOn.{v, u, w} D where
  functor := (Functor.const Dᵒᵖ).obj (AddCommMonCat.of M)
  charInj _ := isCharacteristicallyInjective_id
  fsmIso _ _ := Function.bijective_id

variable {D : Type u} [Category.{v} D]

@[simp] theorem MonoidOn.const_val (M : Type w) [AddCommMonoid M] (A : D) :
    (MonoidOn.const D M).val A = M := rfl

@[simp] theorem MonoidOn.const_map (M : Type w) [AddCommMonoid M] {A B : D} (α : B ⟶ A)
    (x : M) : (MonoidOn.const D M).map α x = x := rfl

/-! ### `𝔽_Φ` が totally epimorphic になる十分条件 -/

/-- `𝒟` の hom が subsingleton で、`Φ(A)` が加法的に簡約的なら、
`𝔽_Φ` は **totally epimorphic** である。

`f ≫ g = f ≫ h` から `deg` は `ℕ+` の簡約性で、`base` は subsingleton で、
`div` は `Φ.map` の単射性(条件 (a))と `Φ(A)` の簡約性で一致する。 -/
theorem isTotallyEpimorphic_elemFrobCat {Φ : MonoidOn.{v, u, w} D}
    (hD : ∀ X Y : D, Subsingleton (X ⟶ Y))
    (hc : ∀ A : D, IsCancelAdd (Φ.val A)) :
    IsTotallyEpimorphic (ElemFrobCat Φ) := by
  intro X Y f
  refine ⟨fun {Z} g h hgh => ?_⟩
  have hdeg : g.deg = h.deg := by
    have := congrArg ElemFrobCat.Hom.deg hgh
    simp only [ElemFrobCat.comp_deg] at this
    exact mul_right_cancel this
  have hdiv : g.div = h.div := by
    have hh := congrArg ElemFrobCat.Hom.div hgh
    simp only [ElemFrobCat.div_comp, hdeg] at hh
    letI := hc Y.base
    exact Φ.map_injective f.base (add_right_cancel hh)
  exact ElemFrobCat.Hom.ext ((hD _ _).elim _ _) hdiv hdeg

/-! ### ★witness 本体 -/

/-- witness の divisor monoid —— `Vee` 上の定数関手 `ℕ`。 -/
abbrev wΦ : MonoidOn.{0, 0, 0} Vee := MonoidOn.const Vee ℕ

/-- witness の `𝒞` —— `𝔽_Φ` 自身。 -/
abbrev wC : Type := ElemFrobCat wΦ

/-- ★**本物の pre-Frobenioid**。`Definition 1.1, (iv)` が課す条件をすべて満たす。 -/
def wP : PreFrobenioid wC wΦ where
  toElem := 𝟭 _
  divisorial _ := isDivisorial_nat
  totEpiC := isTotallyEpimorphic_elemFrobCat (fun X Y => Preorder.subsingleton_hom X Y)
    (fun _ => inferInstanceAs (IsCancelAdd ℕ))
  totEpiD := isTotallyEpimorphic_vee

@[simp] theorem wP_degFr {A B : wC} (φ : A ⟶ B) : wP.degFr φ = φ.deg := rfl
@[simp] theorem wP_Div {A B : wC} (φ : A ⟶ B) : wP.Div φ = φ.div := rfl
@[simp] theorem wP_Base {A B : wC} (φ : A ⟶ B) : wP.Base φ = φ.base := rfl

/-- witness の対象 —— `Vee` の `top` の上。 -/
abbrev wTop : wC := ⟨Vee.top⟩
/-- witness の対象 —— `Vee` の `left` の上。 -/
abbrev wLeft : wC := ⟨Vee.left⟩

/-- `left ⟶ top` の射(`Vee` の中)。 -/
abbrev vLeftTop : Vee.left ⟶ Vee.top := homOfLE (Or.inr rfl)

/-! ### ★非退化 —— `Definition 1.2` の各語に両方向の例

`wTop` の自己射は `(base, div, deg) = (𝟙, d, n)` の形で、`d : ℕ`、`n : ℕ+` が自由に動く。
`linear` は `n`、`isometric` は `d` を見るので、両方向の例がここで取れる。 -/

/-- 次数 2 の自己射(零因子 0)。 -/
def wDeg2 : wTop ⟶ wTop := ⟨𝟙 _, (0 : ℕ), 2⟩
/-- 零因子 1 の自己射(次数 1)。 -/
def wDiv1 : wTop ⟶ wTop := ⟨𝟙 _, (1 : ℕ), 1⟩
/-- `wLeft ⟶ wTop` の射(底が同型でない)。 -/
def wNonBaseIso : wLeft ⟶ wTop := ⟨vLeftTop, (0 : ℕ), 1⟩

/-- 恒等射は **linear**。 -/
theorem wLinear_id : IsLinear wP (𝟙 wTop) := isLinear_id wP wTop

/-- ★次数 2 の射は **linear でない**。 -/
theorem wNot_linear_deg2 : ¬ IsLinear wP wDeg2 := by
  intro h
  have h2 : (2 : ℕ+) = 1 := h
  exact absurd h2 (by decide)

/-- 恒等射は **isometric**。 -/
theorem wIsometric_id : IsIsometric wP (𝟙 wTop) := isIsometric_id wP wTop

/-- ★零因子 1 の射は **isometric でない**。 -/
theorem wNot_isometric_div1 : ¬ IsIsometric wP wDiv1 := by
  intro h
  have h2 : (1 : ℕ) = 0 := h
  exact absurd h2 (by decide)

/-- `Vee` の `left ⟶ top` は同型でない(`top ≤ left` が成り立たない)。 -/
theorem not_isIso_vLeftTop : ¬ IsIso vLeftTop := by
  intro h
  have : Vee.top ≤ Vee.left := leOfHom (inv vLeftTop)
  rcases this with h1 | h1 <;> exact absurd h1 (by decide)

/-- 恒等射は **base-isomorphism**。 -/
theorem wBaseIso_id : IsBaseIsomorphism wP (𝟙 wTop) :=
  isBaseIsomorphism_of_isIso wP (𝟙 wTop)

/-- ★`wLeft ⟶ wTop` は **base-isomorphism でない**。

★ここで底の圏を `Vee` にした意味が出る——一射圏を使っていたら
`base-isomorphism` は恒真になり、この反例が取れなかった。 -/
theorem wNot_baseIso : ¬ IsBaseIsomorphism wP wNonBaseIso := not_isIso_vLeftTop

/-- 恒等射は **pre-step**。 -/
theorem wPreStep_id : IsPreStep wP (𝟙 wTop) := isPreStep_id wP wTop

/-- ★`wLeft ⟶ wTop` は **pre-step でない**(base-isomorphism でないから)。 -/
theorem wNot_preStep : ¬ IsPreStep wP wNonBaseIso := fun h => wNot_baseIso h.2

/-- ★次数 2 の射も **pre-step でない**(linear でないから)。 -/
theorem wNot_preStep_deg2 : ¬ IsPreStep wP wDeg2 := fun h => wNot_linear_deg2 h.1

/-- ★**`step` が空でない** —— 零因子 1 の射は step(同型でない pre-step)。

★これが `Definition 1.2, (iii)` の核心である: `pre-step` のうち
「零因子が 0 でないもの」が本当の `step` になる。 -/
theorem wStep_div1 : IsStep wP wDiv1 := by
  have hbase : IsBaseIsomorphism wP wDiv1 := by
    show IsIso (𝟙 wTop.base)
    infer_instance
  refine ⟨⟨rfl, hbase⟩, fun h => ?_⟩
  have h2 : (1 : ℕ) = 0 :=
    ElemFrobCat.div_eq_zero_of_isIso (fun _ => isSharp_nat) wDiv1
  exact absurd h2 (by decide)

/-- ★**この witness では、すべての対象が `isotropic` である。**

isometric(`div = 0`)な pre-step(`deg = 1`、底が同型)は、成分がすべて
恒等射の値になるので同型である。★`isotropic` の**満たす例**。 -/
theorem wIsotropic (A : wC) : IsIsotropic wP A := by
  intro Dd φ hiso hstep
  have hdeg : ElemFrobCat.Hom.deg φ = 1 := hstep.1
  have hdiv : ElemFrobCat.Hom.div φ = 0 := hiso
  have hbase : IsIso (ElemFrobCat.Hom.base φ) := hstep.2
  refine ⟨⟨inv (ElemFrobCat.Hom.base φ), 0, 1⟩, ?_, ?_⟩
  · refine ElemFrobCat.Hom.ext ?_ ?_ ?_
    · exact IsIso.hom_inv_id _
    · rw [ElemFrobCat.div_comp, hdiv, smul_zero, add_zero]
      rfl
    · simp [hdeg]
  · refine ElemFrobCat.Hom.ext ?_ ?_ ?_
    · exact IsIso.inv_hom_id _
    · rw [ElemFrobCat.div_comp, hdiv, smul_zero, add_zero]
      rfl
    · simp [hdeg]

/-- ★`isotropic` が空虚な条件でないことの確認 —— `wDiv1` は
**isometric でない** pre-step なので、`isotropic` の条件の前提に入らない。
つまり条件は「isometric pre-step 全体」に本当に効いている。 -/
theorem wDiv1_not_in_isotropic_hypothesis : ¬ IsIsometric wP wDiv1 :=
  wNot_isometric_div1

/-! ### ★`Definition 1.1` 側の非退化 —— pre-Frobenioid が空でないこと -/

/-- ★**`𝒞` は totally epimorphic だが `Type` はそうでない**——
`Definition 1.1, (iv)` の条件が自明でないことの確認。 -/
theorem wP_totEpi_nontrivial :
    IsTotallyEpimorphic wC ∧ ¬ IsTotallyEpimorphic (Type) :=
  ⟨wP.totEpiC, not_isTotallyEpimorphic_type⟩

/-- ★**`Φ` が divisorial であることは自明でない**——`ℕ` は divisorial だが
`ℤ` は(sharp でないので)divisorial でない。 -/
theorem wP_divisorial_nontrivial : IsDivisorial ℕ ∧ ¬ IsDivisorial ℤ :=
  ⟨isDivisorial_nat, not_isDivisorial_int⟩

end ABC3.Found.FrdI
