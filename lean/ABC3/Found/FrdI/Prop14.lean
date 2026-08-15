import ABC3.Found.FrdI.NegControl

/-!
# [FrdI] Proposition 1.4 —— co-angular と LB-invertible な射

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、
物理 p.25–p.26(**400 dpi 目視確認 2026-08-15。証明本文まで読むのは今回が初めて**)。

原文 (FrdI p.25):
> Proposition 1.4. (Co-angular and LB-invertible Morphisms) Let Φ be

原文 (FrdI p.25):
> a divisorial monoid on a connected, totally epimorphic category D; C →FΦ a

## ★主張は5つある((i)(ii) だけではない)

原文 (FrdI p.25):
> (i) Suppose that the codomain of any arrow of C whose domain is equal to A is

原文 (FrdI p.25):
> (ii) Suppose that C is a Frobenioid. Then φ is a pull-back morphism if

原文 (FrdI p.26):
> (iii) Suppose that C is a Frobenioid. Then every LB-invertible pre-step

原文 (FrdI p.26):
> (iv) Suppose that C is a Frobenioid. Then a morphism φ of C is co-angular

原文 (FrdI p.26):
> (v) Suppose that C is a Frobenioid. Then a morphism φ of C is LB-invertible

## ★仮定の強さが違う

* **(i)** は `pre-Frobenioid` のまま。追加の仮定は「`A` から出る射の終域がすべて isotropic」
* **(ii)–(v)** は `𝒞` が **Frobenioid** であることを要求する

★したがって (i) だけは `FrobenioidCore` を引数に取らない。

## ★原文の証明の段数(目視で数えた)

原文 (FrdI p.26):
> Proof. Assertion (i) follows formally from the definitions of the terms “isotropic”,

| 主張 | 原文の段 | 引く先 |
|---|---|---|
| (i) | **1段**(「定義から形式的に従う」) | Definition 1.2, (i)(iii)(iv) |
| (ii) → | **1段** | Definition 1.3, (iv), (b) |
| (ii) ← | **2段**(分解を取る / co-angular で `β` を潰す) | Definition 1.3, (iv), (a); Remark 1.1.1 |
| (iii) | **1段**(2つの道が示される) | Definition 1.3, (v)(b) **または** (ii) |
| (iv) → | **1段** | Definition 1.3, (iii), (a) |
| (iv) ← | **3段**(`β` を3分解 / `𝒟` の totally epimorphic / co-angular) | Remark 1.1.1; §0 |
| (v) → | **2段** | Remark 1.1.1; Definition 1.3, (iii), (a) |
| (v) ← | **3段** | (iv); Remark 1.1.1; (iii) |

★**合計 14 段**。★**本ファイルは (i)–(v) の全 14 段を写した。止まった段は無い。**

## ★`Skeleton/` か `Found/` か —— **`Found/` を選んだ**

判断の根拠は**原文の段数**である:

* (i) は原文自身が「follows formally from the definitions」と書いており、**1段**。
  実際 Lean でも1行で済んだ。
* (iii) は `Definition 1.3, (ii)` の本質的一意性に当てるだけで、**1段**。
* (ii) は分解を取って co-angular で潰す **3段**。
* (iv)(v) は最も長いが、各段が `FrobenioidCore` の条項を名指しで引くだけである。

いずれも「証明の見通しが立たない」状態ではないので、`sorry` を置く理由が無い。

★**(iii) の道の選択**: 原文は「`Definition 1.3, (v), (b)` の分解の一意性、**または**
`Definition 1.3, (ii)` の Frobenius 型射の本質的一意性から」と2つの道を示す。
本ファイルは**後者**を採った。★前者の道は、自己監査で足した `preStepFactorUniq`
(一意性節)が無ければ通らなかった——**緩い転写が下流で効かなくなる実例**である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

universe v u w

variable {D : Type u} [Category.{v} D] {C : Type u} [Category.{v} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ### ★(i) —— pre-Frobenioid のままで示せる

原文 (FrdI p.25):
> isotropic. Then φ is co-angular. In particular, φ is a morphism of Frobenius
-/

/-- **[FrdI] Proposition 1.4, (i)** —— `A` から出る射の終域がすべて isotropic なら、
`A` を始域とする射は co-angular。

★証明は原文どおり**1段**: co-angular の条件に現れる `β` は
「`A` から出る射 `γ` の終域 `X`」から出る isometric pre-step なので、
`X` が isotropic であることから直ちに同型。 -/
theorem prop_1_4_i {A B : C} (φ : A ⟶ B)
    (hiso : ∀ (X : C), (A ⟶ X) → IsIsotropic P X) : IsCoAngular P φ :=
  fun X _ γ β _ _ _ hβiso hβstep _ => hiso X γ _ β hβiso hβstep

/-- **[FrdI] Proposition 1.4, (i)** の「In particular」。

原文 (FrdI p.25):
> type if and only if it is an isometric base-isomorphism.

★`Frobenius type = LB-invertible ∧ base-iso = (co-angular ∧ isometric) ∧ base-iso` で、
(i) により co-angular が自動なので、残るのは `isometric ∧ base-iso`。 -/
theorem prop_1_4_i_frobeniusType {A B : C} (φ : A ⟶ B)
    (hiso : ∀ (X : C), (A ⟶ X) → IsIsotropic P X) :
    IsFrobeniusType P φ ↔ (IsIsometric P φ ∧ IsBaseIsomorphism P φ) :=
  ⟨fun h => ⟨h.1.2, h.2⟩, fun h => ⟨⟨prop_1_4_i P φ hiso, h.1⟩, h.2⟩⟩

/-! ### ★(iii) —— Frobenioid が要る

原文 (FrdI p.26):
> is an isomorphism.
-/

/-- **[FrdI] Proposition 1.4, (iii)** —— LB-invertible な pre-step は同型。

★原文は2つの道を示すが、ここでは後者
(`Definition 1.3, (ii)` の Frobenius 型射の本質的一意性)を採る。

`φ` は LB-invertible な pre-step なので**次数 1 の Frobenius 型射**であり、
`𝟙 A` もそうである。本質的一意性から同型 `β : B ⟶ A` が `φ ≫ β = 𝟙 A` を満たし、
`β` が同型なので `φ = inv β` は同型。 -/
theorem prop_1_4_iii (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hlb : IsLBInvertible P φ) (hst : IsPreStep P φ) : IsIso φ := by
  have hft : IsFrobeniusType P φ := ⟨hlb, hst.2⟩
  have hid : IsFrobeniusType P (𝟙 A) :=
    ⟨⟨isCoAngular_id P A, isIsometric_id P A⟩, isBaseIsomorphism_of_isIso P (𝟙 A)⟩
  have hdeg : P.degFr φ = P.degFr (𝟙 A) := by
    rw [P.degFr_id]
    exact hst.1
  obtain ⟨β, hβiso, hβ⟩ := F.frobDegUniq A B A φ (𝟙 A) hft hid hdeg
  haveI := hβiso
  have hφ : φ = 𝟙 A ≫ inv β := (IsIso.eq_comp_inv β).mpr hβ
  rw [hφ]
  infer_instance

/-! ### ★(ii) —— Frobenioid が要る

原文 (FrdI p.25):
> and only if it is an LB-invertible linear morphism [i.e., a co-angular linear

原文 (FrdI p.26):
> that β is an isomorphism, hence that φ is a pull-back morphism, as desired. As-
-/

/-- **(ii) の順方向** —— pull-back morphism なら LB-invertible かつ linear。

★これは `Definition 1.3, (iv), (b)` そのもの(`FrobenioidCore.pullBackLB`)。**1段**。 -/
theorem prop_1_4_ii_mp (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (h : IsPullBack P φ) : IsLBInvertible P φ ∧ IsLinear P φ := F.pullBackLB φ h

/-- **(ii) の逆方向** —— LB-invertible かつ linear なら pull-back morphism。

★原文の3段をそのまま写す:
1. `Definition 1.3, (iv), (a)` の分解 `φ = γ ≫ β ≫ α` を取る
2. `φ` が linear かつ isometric であることから、`γ` は同型、`γ ≫ β` は isometric pre-step
3. `φ` が co-angular なので `γ ≫ β` は同型、よって `φ` は pull-back morphism -/
theorem prop_1_4_ii_mpr (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hlb : IsLBInvertible P φ) (hlin : IsLinear P φ) : IsPullBack P φ := by
  obtain ⟨X, Y, γ, β, α, hfac, hγ, hβ, hα⟩ := F.arbFactor φ
  obtain ⟨hαlb, hαlin⟩ := F.pullBackLB α hα
  -- 段2: `γ` の次数は 1
  have hdγ : P.degFr γ = 1 := by
    have h1 : P.degFr φ = 1 := hlin
    rw [hfac, P.degFr_comp, P.degFr_comp, show P.degFr α = 1 from hαlin,
      show P.degFr β = 1 from hβ.1, one_mul, one_mul] at h1
    exact h1
  -- 次数 1 の Frobenius 型射は pre-step でもあるので、(iii) により同型
  haveI hγiso : IsIso γ := prop_1_4_iii P F γ hγ.1 ⟨hdγ, hγ.2⟩
  -- 段2 続き: `γ ≫ β` は pre-step で、しかも isometric
  have hstep : IsPreStep P (γ ≫ β) := by
    refine ⟨?_, isBaseIsomorphism_comp P (isBaseIsomorphism_of_isIso P γ) hβ.2⟩
    show P.degFr (γ ≫ β) = 1
    rw [P.degFr_comp, show P.degFr β = 1 from hβ.1, hdγ, one_mul]
  have hisom : IsIsometric P (γ ≫ β) := by
    have h0 : P.Div φ = 0 := hlb.2
    rw [hfac, ← Category.assoc, P.Div_comp, show P.Div α = 0 from hαlb.2,
      show P.degFr α = 1 from hαlin, map_zero, zero_add] at h0
    show P.Div (γ ≫ β) = 0
    simpa using h0
  -- 段3: co-angular を分解 `φ = 𝟙 ≫ (γ ≫ β) ≫ α` に当てる
  haveI : IsIso (γ ≫ β) := by
    refine hlb.1 A Y (𝟙 A) (γ ≫ β) α ?_ hαlin hisom hstep
      (Or.inr (isBaseIsomorphism_of_isIso P (𝟙 A)))
    rw [Category.id_comp]
    exact hfac.trans (Category.assoc _ _ _).symm
  -- `φ = (γ ≫ β) ≫ α` で、前半は同型、後半は pull-back
  have hfac2 : φ = (γ ≫ β) ≫ α := hfac.trans (Category.assoc _ _ _).symm
  rw [hfac2]
  exact IsPullBack.comp P (isPullBack_of_isIso P (γ ≫ β)) hα

/-- **[FrdI] Proposition 1.4, (ii)** —— 両方向。 -/
theorem prop_1_4_ii (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B) :
    IsPullBack P φ ↔ (IsLBInvertible P φ ∧ IsLinear P φ) :=
  ⟨prop_1_4_ii_mp P F φ, fun h => prop_1_4_ii_mpr P F φ h.1 h.2⟩

/-! ### ★(iv) —— Frobenioid が要る

原文 (FrdI p.26):
> if and only if, in the factorization φ = α ◦β ◦γ of Definition 1.3, (iv), (a), the

原文 (FrdI p.26):
> Next, we consider assertion (iv). If β is co-angular, then since α, γ are co-
-/

/-- **(iv) の順方向** —— 分解の `β` が co-angular なら `φ` も co-angular。

★原文どおり**1段**: `α` は pull-back なので (ii) により co-angular、
`γ` は Frobenius 型なので定義から co-angular。あとは
`Definition 1.3, (iii), (a)`(co-angular は合成で閉じる)を2回当てるだけ。 -/
theorem prop_1_4_iv_mp (F : FrobenioidCore P) {A B X Y : C} {φ : A ⟶ B}
    {γ : A ⟶ X} {β : X ⟶ Y} {α : Y ⟶ B} (hfac : φ = γ ≫ β ≫ α)
    (hγ : IsFrobeniusType P γ) (hα : IsPullBack P α) (hβco : IsCoAngular P β) :
    IsCoAngular P φ := by
  have hαco : IsCoAngular P α := (F.pullBackLB α hα).1.1
  rw [hfac]
  exact F.coAngularComp γ (β ≫ α) hγ.1.1 (F.coAngularComp β α hβco hαco)

/-- **(iv) の逆方向** —— `φ` が co-angular なら分解の `β` も co-angular。

★原文の3段をそのまま写す:
1. `β` の3分解 `β = γ' ≫ β' ≫ α'` に対し、次数の乗法性から `γ'`、`α'` も linear
2. `𝒟` が totally epimorphic なので、`Base β` の同型性から `Base γ'`、`Base α'` も同型
   (`isIso_of_isIso_comp` = 原文 p.16 の「in passing」の観察)
3. `φ = (γ ≫ γ') ≫ β' ≫ (α' ≫ α)` に `φ` の co-angular 性を当てて `β'` が同型 -/
theorem prop_1_4_iv_mpr (F : FrobenioidCore P) {A B X Y : C} {φ : A ⟶ B}
    {γ : A ⟶ X} {β : X ⟶ Y} {α : Y ⟶ B} (hfac : φ = γ ≫ β ≫ α)
    (hγ : IsFrobeniusType P γ) (hβ : IsPreStep P β) (hα : IsPullBack P α)
    (hφco : IsCoAngular P φ) : IsCoAngular P β := by
  intro X' Y' γ' β' α' hfac' hα'lin hβ'iso hβ'step _
  obtain ⟨hαlb, hαlin⟩ := F.pullBackLB α hα
  -- 段1: 次数の乗法性から `γ'` も linear
  have hdβ : P.degFr β = 1 := hβ.1
  have hdγ' : P.degFr γ' = 1 := by
    rw [hfac', P.degFr_comp, P.degFr_comp, show P.degFr α' = 1 from hα'lin,
      show P.degFr β' = 1 from hβ'step.1, one_mul, one_mul] at hdβ
    exact hdβ
  -- 段2: `𝒟` の totally epimorphic 性から `Base γ'`、`Base α'` が同型
  haveI hbβ : IsIso (P.Base β) := hβ.2
  haveI : IsIso (P.Base γ' ≫ P.Base β' ≫ P.Base α') := by
    rw [← P.Base_comp, ← P.Base_comp, ← hfac']
    exact hbβ
  have h2 := isIso_of_isIso_comp P.totEpiD (P.Base γ') (P.Base β' ≫ P.Base α')
  haveI : IsIso (P.Base β' ≫ P.Base α') := h2.2
  have h3 := isIso_of_isIso_comp P.totEpiD (P.Base β') (P.Base α')
  -- `γ'` と `α'` は pre-step
  have hγ'step : IsPreStep P γ' := ⟨hdγ', h2.1⟩
  have hα'step : IsPreStep P α' := ⟨hα'lin, h3.2⟩
  -- 段3: `φ = (γ ≫ γ') ≫ β' ≫ (α' ≫ α)` に co-angular を当てる
  refine hφco X' Y' (γ ≫ γ') β' (α' ≫ α) ?_ ?_ hβ'iso hβ'step (Or.inr ?_)
  · rw [hfac, hfac']
    simp
  · show P.degFr (α' ≫ α) = 1
    rw [P.degFr_comp, show P.degFr α = 1 from hαlin, show P.degFr α' = 1 from hα'lin, one_mul]
  · exact isBaseIsomorphism_comp P hγ.2 hγ'step.2

/-! ### ★(v) —— Frobenioid が要る

原文 (FrdI p.26):
> if and only if it is of the form α ◦β, where α is a pull-back morphism, and β is

原文 (FrdI p.26):
> Finally, we consider assertion (v). If φ = α◦β, where α is a pull-back morphism,
-/

/-- **(v) の順方向** —— `pull-back ∘ Frobenius 型` は LB-invertible。

★原文どおり**2段**: `α`、`β` がともに LB-invertible なので、
Remark 1.1.1 から `φ` は isometric、`Definition 1.3, (iii), (a)` から co-angular。 -/
theorem prop_1_4_v_mp (F : FrobenioidCore P) {A X B : C} {φ : A ⟶ B}
    {β : A ⟶ X} {α : X ⟶ B} (hfac : φ = β ≫ α)
    (hβ : IsFrobeniusType P β) (hα : IsPullBack P α) : IsLBInvertible P φ := by
  obtain ⟨hαlb, hαlin⟩ := F.pullBackLB α hα
  constructor
  · rw [hfac]
    exact F.coAngularComp β α hβ.1.1 hαlb.1
  · show P.Div φ = 0
    rw [hfac, P.Div_comp, show P.Div α = 0 from hαlb.2,
      show P.Div β = 0 from hβ.1.2, map_zero, zero_add, smul_zero]

/-- **(v) の逆方向** —— LB-invertible なら `pull-back ∘ Frobenius 型` の形。

★原文の3段をそのまま写す:
1. `Definition 1.3, (iv), (a)` の分解を取り、(iv) により `β` は co-angular
2. Remark 1.1.1 により `β` は isometric(`Φ.map` の単射性を使う)
3. よって `β` は LB-invertible な pre-step、(iii) により同型。
   `γ ≫ β` は Frobenius 型なので `φ = (γ ≫ β) ≫ α` が求める形 -/
theorem prop_1_4_v_mpr (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hlb : IsLBInvertible P φ) :
    ∃ (X : C) (β : A ⟶ X) (α : X ⟶ B),
      φ = β ≫ α ∧ IsFrobeniusType P β ∧ IsPullBack P α := by
  obtain ⟨X, Y, γ, β, α, hfac, hγ, hβ, hα⟩ := F.arbFactor φ
  obtain ⟨hαlb, hαlin⟩ := F.pullBackLB α hα
  -- 段1: (iv) により `β` は co-angular
  have hβco : IsCoAngular P β := prop_1_4_iv_mpr P F hfac hγ hβ hα hlb.1
  -- 段2: `β` は isometric
  have hβiso : IsIsometric P β := by
    have h0 : P.Div φ = 0 := hlb.2
    rw [hfac, P.Div_comp, show P.Div γ = 0 from hγ.1.2, smul_zero, add_zero] at h0
    have h1 : P.Div (β ≫ α) = 0 := Φ.map_injective (P.Base γ) (by simpa using h0)
    rw [P.Div_comp, show P.Div α = 0 from hαlb.2, show P.degFr α = 1 from hαlin,
      map_zero, zero_add] at h1
    show P.Div β = 0
    simpa using h1
  -- 段3: `β` は LB-invertible な pre-step、よって同型
  haveI : IsIso β := prop_1_4_iii P F β ⟨hβco, hβiso⟩ hβ
  refine ⟨Y, γ ≫ β, α, hfac.trans (Category.assoc _ _ _).symm, ?_, hα⟩
  refine ⟨⟨F.coAngularComp γ β hγ.1.1 hβco, ?_⟩,
    isBaseIsomorphism_comp P hγ.2 hβ.2⟩
  show P.Div (γ ≫ β) = 0
  rw [P.Div_comp, show P.Div β = 0 from hβiso, show P.Div γ = 0 from hγ.1.2,
    map_zero, zero_add, smul_zero]

/-! ### ★模型 `wP` の上での確認

★**模型で成り立つことは一般の証明ではない。** ここでやるのは
「一般定理を模型に当てた結果が、模型を**独立に**調べて得た特徴づけと一致するか」
という**転写の検査**である。 -/

/-- (ii) を模型 `wP` に当てたもの。 -/
theorem wProp_1_4_ii {A B : wC} (φ : A ⟶ B) :
    IsPullBack wP φ ↔ (IsLBInvertible wP φ ∧ IsLinear wP φ) :=
  prop_1_4_ii wP wIsFrobenioid.core φ

/-- ★**転写の検査** —— 一般定理 (ii) の右辺が、模型を独立に調べて得た
`wIsPullBack_iff`(`wd φ = 0 ∧ wn φ = 1`)と**一致する**。

★2つは別々に得たものである: 左辺は `Definition 1.3` の条項から一般に証明したもの、
右辺は `𝔽_Φ` の合成則を直接計算して得たもの。**一致は転写が正しいことの証拠**である。 -/
theorem wProp_1_4_ii_agrees {A B : wC} (φ : A ⟶ B) :
    (IsLBInvertible wP φ ∧ IsLinear wP φ) ↔ (wd φ = 0 ∧ wn φ = 1) :=
  (wProp_1_4_ii φ).symm.trans (wIsPullBack_iff φ)

/-- (i) を模型 `wP` に当てたもの —— `wP` は全対象 isotropic なので仮定は自動で満たされる。 -/
theorem wProp_1_4_i {A B : wC} (φ : A ⟶ B) : IsCoAngular wP φ :=
  prop_1_4_i wP φ (fun X _ => wIsotropic X)

/-! ### ★★負の対照 —— (i) の仮定は効いている

★`Proposition 1.4, (i)` は「`A` から出る射の終域がすべて isotropic」を仮定する。
**その仮定を外すと結論(co-angular)が落ちる**ことを、実物で示す。

★`𝔽_Φ` とその広い部分圏では**どの対象も isotropic** である
(isometric pre-step は `(b, 0, 1)` で `b` が同型、したがって既に同型)。
そこで **`𝒞 → 𝔽_Φ` を定数関手に取る**——すると `Base` はつねに `𝟙`、
`Div` はつねに `0`、`deg_Fr` はつねに `1` になり、
**`𝒞` のすべての射が isometric pre-step になる**。 -/

/-- 負の対照の `𝒞 → 𝔽_Φ` —— 定数関手。 -/
def cToElem : Vee ⥤ wC := (Functor.const Vee).obj wTop

/-- ★負の対照の pre-Frobenioid —— `𝒞 = Vee`、`𝒞 → 𝔽_Φ` は定数関手。 -/
def cP : PreFrobenioid Vee wΦ where
  toElem := cToElem
  divisorial _ := isDivisorial_nat
  totEpiC := isTotallyEpimorphic_vee
  totEpiD := isTotallyEpimorphic_vee

@[simp] theorem cP_Base {A B : Vee} (φ : A ⟶ B) : cP.Base φ = 𝟙 Vee.top := rfl
@[simp] theorem cP_Div {A B : Vee} (φ : A ⟶ B) : cP.Div φ = 0 := rfl
@[simp] theorem cP_degFr {A B : Vee} (φ : A ⟶ B) : cP.degFr φ = 1 := rfl

/-- 定数関手のもとでは、`Vee` の**すべての射**が isometric pre-step。 -/
theorem cP_isometric {A B : Vee} (φ : A ⟶ B) : IsIsometric cP φ := rfl

theorem cP_preStep {A B : Vee} (φ : A ⟶ B) : IsPreStep cP φ :=
  ⟨rfl, (inferInstance : IsIso (𝟙 Vee.top))⟩

/-- ★**`Vee.left` は isotropic でない** —— `left ⟶ top` は isometric pre-step だが同型でない。 -/
theorem cP_not_isotropic_left : ¬ IsIsotropic cP Vee.left := by
  intro h
  exact not_isIso_vLeftTop (h Vee.top vLeftTop (cP_isometric _) (cP_preStep _))

/-- ★★**(i) の仮定を外すと結論が落ちる** —— `left ⟶ top` は co-angular でない。

★したがって `Proposition 1.4, (i)` の isotropy 仮定は**効いている**
(load-bearing である)。 -/
theorem cP_not_coAngular : ¬ IsCoAngular cP vLeftTop := by
  intro h
  exact not_isIso_vLeftTop
    (h Vee.left Vee.top (𝟙 Vee.left) vLeftTop (𝟙 Vee.top) (by simp) rfl
      (cP_isometric _) (cP_preStep _)
      (Or.inr (isBaseIsomorphism_of_isIso cP (𝟙 Vee.left))))

/-- ★**負の対照のまとめ** —— (i) は仮定つきでのみ成り立つ。

満たす側: `wP` では全対象が isotropic で、全射が co-angular。
満たさない側: `cP` では `left` が isotropic でなく、`left ⟶ top` が co-angular でない。 -/
theorem prop_1_4_i_hypothesis_is_load_bearing :
    (∀ {A B : wC} (φ : A ⟶ B), IsCoAngular wP φ) ∧ ¬ IsCoAngular cP vLeftTop :=
  ⟨fun φ => wProp_1_4_i φ, cP_not_coAngular⟩

/-! ### ★出典の紐付け(`.src`) -/

def prop_1_4_i.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 25, item := "Proposition 1.4, (i)",
    sectionId := "frdi-prop-1-4-i" }

def prop_1_4_ii.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 25, item := "Proposition 1.4, (ii)",
    sectionId := "frdi-prop-1-4-ii" }

def prop_1_4_iii.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 26, item := "Proposition 1.4, (iii)",
    sectionId := "frdi-prop-1-4-iii" }

def prop_1_4_iv_mp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 26, item := "Proposition 1.4, (iv)",
    sectionId := "frdi-prop-1-4-iv" }

def prop_1_4_v_mp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 26, item := "Proposition 1.4, (v)",
    sectionId := "frdi-prop-1-4-v" }

def cP.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 25, item := "Proposition 1.4, (i) — 負の対照",
    sectionId := "frdi-prop-1-4-i" }

end ABC3.Found.FrdI
