/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.TrivSecNorm
import ABC3.Found.Arakelov.ArcLocalMetric
import ABC3.Meta.Claim

/-!
# **自明化の側で書いた計量**とそのテンソル積（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★台帳 `arakelov-coherent-base-metric` の段 3

| 段 | 内容 | 場所 |
|---|---|---|
| 1 | `trivValue` がテンソル積で掛け算になる | `TrivTensor.lean` |
| 2 | 切断のノルムが掛け算になる（ファイバー迂回） | `TrivSecNorm.lean` |
| 3 | **基準ノルムを担ぐ構造体と、自明化に依らないこと** | ★★**本ファイル** |

## ★★★★★★★計量を「基準ノルム」で書く

局所自明な前層加群の計量は、自明化ごとの**基準ノルム**

    `h_{V,e}(p) ≝ ‖e⁻¹(1)‖(p) > 0`

で表せる。★そこでは切断のノルムが

    `‖s‖(p) = ‖evalOn p V hp (trivValue F V e s)‖ · h_{V,e}(p)`

と書け、★★**自明化を取り替えると `trivValue` は単元 `u` 倍・`h` は `|u|` 分の 1 倍**
になるので、積は変わらない。★★★それが `normOf_indep` である。

## ★★★★★★★★★そしてテンソル積で掛け算になる

`h_{A⊗B, eA⊗eB} = h_{A,eA} · h_{B,eB}` を満たす計量（`IsTensorOf`）に対して

    `‖s ⊗ t‖ = ‖s‖ · ‖t‖`

★これが **`Classical.choice` の基準計量では落ちていた加法性**である。
★★`normOf_indep` と併せると、**テンソル自明化の上での値だけで計量が決まる**。

## ★残っている段（明示）

★★**存在**——`IsTensorOf` を満たす `LocalMetric X (A ⊗ B)` を `mA`, `mB` から
実際に作る段。★そこでは「`(A⊗B)|_V` は自明だが `A|_V` は自明でない `V`」が問題になり、
両方が自明になる小さい `W ∋ p` へ降りて貼り合わせる必要がある。
★★★値が選択に依らないことは本ファイルの `normOf_indep` と `compat` が与える。

★★★★もう 1 つ: 制限との両立（`h` を `W ⊆ V` へ落とす欄）を構造体に足す段。
本ファイルは**同じ `V` の上での自明化の取り替え**だけを扱う。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite
open ABC3.Found.GenEll

variable {X : Scheme.{0}}

/-! ## ★★自明化の**遷移単元**（正準な形）

★★★`trivValue_congr'` は `∃ u, …` を与えるだけなので、
それを構造体の欄に使うと**「同じ関係を満たす `u`, `u'` が違いうる」**という
空虚の危険が残る（`trivValue F V e` は `Γ(X,V)` の上へ**全射とは限らない**——
`secOn` は大域切断からの制限だから）。
★したがって遷移単元を**式で**定める。 -/

/-- ★★★**自明化 `e` から `e'` への遷移単元** `u = θ 1`。 -/
noncomputable def transUnit (F : X.PresheafOfModules) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V)) :
    ((X.presheaf.obj (op V)) : Type) :=
  ((trivEquiv F V e).symm.trans (trivEquiv F V e')) 1

theorem isUnit_transUnit (F : X.PresheafOfModules) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V)) :
    IsUnit (transUnit F V e e') :=
  linearEquiv_self_isUnit _

/-- ★★**遷移単元は実際に `trivValue` を移す**（`trivValue_congr'` の正準版）。 -/
theorem trivValue_transUnit (F : X.PresheafOfModules) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (s : F.obj (op ⊤)) :
    trivValue F V e' s = trivValue F V e s * transUnit F V e e' := by
  set θ : ((Γ(X, V)) : Type) ≃ₗ[((Γ(X, V)) : Type)] ((Γ(X, V)) : Type) :=
    (trivEquiv F V e).symm.trans (trivEquiv F V e') with hθ
  have h := linearEquiv_self_apply θ (trivEquiv F V e (secOn F V s))
  have h2 : θ (trivEquiv F V e (secOn F V s)) = trivEquiv F V e' (secOn F V s) := by
    rw [hθ]; simp
  rw [trivValue_eq_trivEquiv, trivValue_eq_trivEquiv, ← h2, h]
  rfl

/-- ★★★★**遷移単元はコサイクルである**: `u_{e→e''} = u_{e→e'} · u_{e'→e''}`。

★機構は `linearEquiv_self_apply`（`θ x = x · θ 1`）だけである。 -/
theorem transUnit_trans (F : X.PresheafOfModules) (V : X.Opens)
    (e e' e'' : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V)) :
    transUnit F V e e'' = transUnit F V e e' * transUnit F V e' e'' := by
  have h1 : ((trivEquiv F V e).symm.trans (trivEquiv F V e'')) 1
      = ((trivEquiv F V e').symm.trans (trivEquiv F V e''))
          (((trivEquiv F V e).symm.trans (trivEquiv F V e')) 1) := by
    simp
  show ((trivEquiv F V e).symm.trans (trivEquiv F V e'')) 1 = _
  rw [h1, linearEquiv_self_apply]
  rfl

/-- ★遷移単元の自明な場合。 -/
theorem transUnit_self (F : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V)) :
    transUnit F V e e = 1 := by
  show ((trivEquiv F V e).symm.trans (trivEquiv F V e)) 1 = 1
  simp

/-- ★単元の値は `0` でない——`evalOn` は環準同型だから。 -/
theorem evalOn_ne_zero_of_isUnit (p : Spec (CommRingCat.of ℂ) ⟶ X) (V : X.Opens)
    (hp : p ⁻¹ᵁ V = ⊤) {u : ((X.presheaf.obj (op V)) : Type)} (hu : IsUnit u) :
    evalOn p V hp u ≠ 0 := by
  obtain ⟨v, hv⟩ := hu.exists_right_inv
  intro h0
  have : (1 : (CommRingCat.of ℂ : Type)) = 0 := by
    rw [← evalOn_one V p hp, ← hv, evalOn_mul, h0, zero_mul]
  exact one_ne_zero this

/-- ★`trivSecNorm` の、遷移単元を受け取る形。

★★これが `LocalMetric.compat` と噛み合う。 -/
theorem trivSecNorm_of_congr (F : X.PresheafOfModules) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) (s : F.obj (op ⊤)) :
    trivSecNorm F V e' p hp s
      = trivSecNorm F V e p hp s * ‖evalOn p V hp (transUnit F V e e')‖ := by
  show ‖evalOn p V hp (trivValue F V e' s)‖ = _
  rw [trivValue_transUnit F V e e' s, evalOn_mul, norm_mul]
  rfl

/-! ## ★★★★★★基準ノルムを担ぐ構造体 -/

/-- ★★★★★★**自明化の側で書いた計量**。

`h V e p` は「自明化 `e` の生成切断 `e⁻¹(1)` の `p` でのノルム」である。

★★`compat` が「自明化を `u` 倍すると `h` は `|u|` で割られる」を言う
——これが切断のノルムを自明化に依らなくする欄である。 -/
structure LocalMetric (X : Scheme.{0}) (F : X.PresheafOfModules) where
  /-- 自明化ごとの基準ノルム。 -/
  h : ∀ (V : X.Opens), ((restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V)) →
    (Spec (CommRingCat.of ℂ) ⟶ X) → ℝ
  /-- ★基準ノルムは正である（生成切断は消えない）。 -/
  pos : ∀ V e p, 0 < h V e p
  /-- ★★**自明化を遷移単元 `u` 倍すると `h` は `|u|` で割られる**。

  ★遷移単元は `transUnit`（式で定めた正準なもの）である
  ——`∃ u, …` の `u` を欄に取ると空虚の危険が残るからである。 -/
  compat : ∀ (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤),
    h V e' p * ‖evalOn p V hp (transUnit F V e e')‖ = h V e p

namespace LocalMetric

variable {F : X.PresheafOfModules}

/-- ★★★**切断のノルム**。 -/
noncomputable def normOf (m : LocalMetric X F) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) (s : F.obj (op ⊤)) : ℝ :=
  trivSecNorm F V e p hp s * m.h V e p

theorem normOf_nonneg (m : LocalMetric X F) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) (s : F.obj (op ⊤)) :
    0 ≤ m.normOf V e p hp s :=
  mul_nonneg (norm_nonneg _) (m.pos V e p).le

/-- ★★★★★★**切断のノルムは自明化の取り方に依らない**。

★機構は「`trivValue` は `u` 倍・`h` は `|u|` 分の 1 倍」——積で消える。 -/
theorem normOf_indep (m : LocalMetric X F) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) (s : F.obj (op ⊤)) :
    m.normOf V e' p hp s = m.normOf V e p hp s := by
  show trivSecNorm F V e' p hp s * m.h V e' p = trivSecNorm F V e p hp s * m.h V e p
  rw [trivSecNorm_of_congr F V e e' p hp s, ← m.compat V e e' p hp]
  ring

/-- ★★ノルムが 0 になるのは `trivValue` が `p` で消えるときだけ。 -/
theorem normOf_eq_zero_iff (m : LocalMetric X F) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) (s : F.obj (op ⊤)) :
    m.normOf V e p hp s = 0 ↔ evalOn p V hp (trivValue F V e s) = 0 := by
  show trivSecNorm F V e p hp s * m.h V e p = 0 ↔ _
  rw [mul_eq_zero]
  constructor
  · rintro (h1 | h2)
    · exact norm_eq_zero.1 h1
    · exact absurd h2 (m.pos V e p).ne'
  · intro h1
    exact Or.inl (by show ‖evalOn p V hp (trivValue F V e s)‖ = 0; rw [h1, norm_zero])

end LocalMetric

/-! ## ★★★★★★★★★テンソル積 -/

/-- ★★★★★★**`m` が `mA`, `mB` のテンソル積であること**。

★基準ノルムがテンソル自明化の上で掛け算になる、というだけの条件である。 -/
def IsTensorOf {A B : X.PresheafOfModules} (mA : LocalMetric X A) (mB : LocalMetric X B)
    (m : LocalMetric X (A ⊗ B)) : Prop :=
  ∀ (V : X.Opens)
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X),
    m.h V (tensorTriv eA eB) p = mA.h V eA p * mB.h V eB p

/-- ★★★★★★★★★**テンソル積の計量では切断のノルムが掛け算になる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★これが `Classical.choice` の基準計量では落ちていた加法性の**本体**である。

★★`normOf_indep` と併せると、**テンソル自明化の上での値だけで計量が決まる**。 -/
theorem normOf_tensor {A B : X.PresheafOfModules} {mA : LocalMetric X A} {mB : LocalMetric X B}
    {m : LocalMetric X (A ⊗ B)} (hm : IsTensorOf mA mB m) (V : X.Opens)
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤)
    (s : A.obj (op ⊤)) (t : B.obj (op ⊤)) :
    m.normOf V (tensorTriv eA eB) p hp (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t)
      = mA.normOf V eA p hp s * mB.normOf V eB p hp t := by
  show trivSecNorm (A ⊗ B) V (tensorTriv eA eB) p hp _ * m.h V (tensorTriv eA eB) p = _
  rw [trivSecNorm_tensor, hm V eA eB p]
  show _ = (trivSecNorm A V eA p hp s * mA.h V eA p) * (trivSecNorm B V eB p hp t * mB.h V eB p)
  ring

/-- ★★★★★★★**任意の自明化で見ても掛け算になる**。

★`normOf_indep` を噛ませただけだが、**これが「計量のテンソル積」の意味**である
——テンソル自明化を選ばずに言える。 -/
theorem normOf_tensor' {A B : X.PresheafOfModules} {mA : LocalMetric X A} {mB : LocalMetric X B}
    {m : LocalMetric X (A ⊗ B)} (hm : IsTensorOf mA mB m) (V : X.Opens)
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V))
    (f : (restrictPresheafFunctor X V).obj (A ⊗ B) ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤)
    (s : A.obj (op ⊤)) (t : B.obj (op ⊤)) :
    m.normOf V f p hp (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t)
      = mA.normOf V eA p hp s * mB.normOf V eB p hp t :=
  (m.normOf_indep V (tensorTriv eA eB) f p hp _).trans (normOf_tensor hm V eA eB p hp s t)

/-- ★★**対数は足し算になる**——Green 関数の側の形。 -/
theorem log_normOf_tensor {A B : X.PresheafOfModules} {mA : LocalMetric X A}
    {mB : LocalMetric X B} {m : LocalMetric X (A ⊗ B)} (hm : IsTensorOf mA mB m) (V : X.Opens)
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤)
    (s : A.obj (op ⊤)) (t : B.obj (op ⊤))
    (hs : mA.normOf V eA p hp s ≠ 0) (ht : mB.normOf V eB p hp t ≠ 0) :
    Real.log (m.normOf V (tensorTriv eA eB) p hp (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t))
      = Real.log (mA.normOf V eA p hp s) + Real.log (mB.normOf V eB p hp t) := by
  rw [normOf_tensor hm V eA eB p hp s t, Real.log_mul hs ht]

/-! ## ★★★★★★★★空虚でないこと —— `Ō_X` の標準計量

★★構造層 `𝟙_` には**基準の自明化**（`restrictPresheafUnit`）があるので、
そこで `h = 1` と置き、他の自明化へは遷移単元で運べばよい。
★★★整合条件はコサイクル（`transUnit_trans`）から出る。 -/

section StructMetric

open scoped Classical

/-- ★構造層の**基準の自明化**。 -/
noncomputable def baseTriv (X : Scheme.{0}) (V : X.Opens) :
    (restrictPresheafFunctor X V).obj (𝟙_ X.PresheafOfModules) ≅ 𝟙_ (PresheafModulesOn X V) :=
  (restrictPresheafUnit (X := X) (U := V)).symm

/-- ★★★★★★★★**構造層 `Ō_X` の標準計量**——`LocalMetric` が空虚でないことの目撃者。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★基準の自明化で `h = 1`、他へは遷移単元の絶対値の逆数で運ぶ。
★★整合条件は `transUnit_trans`（コサイクル）そのものである。 -/
noncomputable def structLocalMetric (X : Scheme.{0}) :
    LocalMetric X (𝟙_ X.PresheafOfModules) where
  h V e p :=
    if hp : p ⁻¹ᵁ V = ⊤ then
      ‖evalOn p V hp (transUnit (𝟙_ X.PresheafOfModules) V (baseTriv X V) e)‖⁻¹
    else 1
  pos V e p := by
    by_cases hp : p ⁻¹ᵁ V = ⊤
    · show 0 < (if hp : p ⁻¹ᵁ V = ⊤ then _ else (1 : ℝ))
      rw [dif_pos hp]
      exact inv_pos.2 (lt_of_le_of_ne (norm_nonneg _)
        (Ne.symm (norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit p V hp
          (isUnit_transUnit _ V (baseTriv X V) e)))))
    · show 0 < (if hp : p ⁻¹ᵁ V = ⊤ then _ else (1 : ℝ))
      rw [dif_neg hp]
      exact one_pos
  compat V e e' p hp := by
    show (if hp : p ⁻¹ᵁ V = ⊤ then
        ‖evalOn p V hp (transUnit (𝟙_ X.PresheafOfModules) V (baseTriv X V) e')‖⁻¹
      else 1) * _ = _
    rw [dif_pos hp]
    show ‖evalOn p V hp (transUnit (𝟙_ X.PresheafOfModules) V (baseTriv X V) e')‖⁻¹
        * ‖evalOn p V hp (transUnit (𝟙_ X.PresheafOfModules) V e e')‖
      = (if hp : p ⁻¹ᵁ V = ⊤ then
          ‖evalOn p V hp (transUnit (𝟙_ X.PresheafOfModules) V (baseTriv X V) e)‖⁻¹
        else 1)
    rw [dif_pos hp, transUnit_trans (𝟙_ X.PresheafOfModules) V (baseTriv X V) e e',
      evalOn_mul, norm_mul, mul_inv]
    have hne : ‖evalOn p V hp (transUnit (𝟙_ X.PresheafOfModules) V e e')‖ ≠ 0 :=
      norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit p V hp (isUnit_transUnit _ V e e'))
    field_simp

/-- ★基準の自明化では `h = 1` である。 -/
theorem structLocalMetric_h_baseTriv (X : Scheme.{0}) (V : X.Opens)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) :
    (structLocalMetric X).h V (baseTriv X V) p = 1 := by
  show (if hp : p ⁻¹ᵁ V = ⊤ then
      ‖evalOn p V hp (transUnit (𝟙_ X.PresheafOfModules) V (baseTriv X V) (baseTriv X V))‖⁻¹
    else 1) = 1
  rw [dif_pos hp, transUnit_self, evalOn_one, norm_one, inv_one]

/-- ★★**したがって `Ō_X` では切断のノルムは `|s(p)|` そのものである**（基準の自明化で）。 -/
theorem structLocalMetric_normOf_baseTriv (X : Scheme.{0}) (V : X.Opens)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤)
    (s : (𝟙_ X.PresheafOfModules).obj (op ⊤)) :
    (structLocalMetric X).normOf V (baseTriv X V) p hp s
      = ‖evalOn p V hp (trivValue (𝟙_ X.PresheafOfModules) V (baseTriv X V) s)‖ := by
  show trivSecNorm _ V (baseTriv X V) p hp s * (structLocalMetric X).h V (baseTriv X V) p = _
  rw [structLocalMetric_h_baseTriv X V p hp, mul_one]
  rfl

end StructMetric

/-! ### ★出典の紐付け(`.src`) -/

def LocalMetric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(自明化の側で書いた計量——基準ノルム h_{V,e})",
    sectionId := "genell-def-1-1-i" }

def LocalMetric.normOf_indep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(切断のノルムが自明化の取り方に依らないこと)",
    sectionId := "genell-def-1-1-i" }

def transUnit_trans.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(自明化の遷移単元はコサイクルであること)",
    sectionId := "genell-def-1-1-i" }

def structLocalMetric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(構造層 Ō_X の標準計量——LocalMetric が空虚でないこと)",
    sectionId := "genell-def-1-1-i" }

def normOf_tensor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(テンソル積の計量では切断のノルムが掛け算になること)",
    sectionId := "genell-def-1-1-i" }

def normOf_tensor.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "trivSecNorm_tensor(切断のノルムがテンソル積で掛け算になる)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.trivSecNorm_tensor") 3,
    .citation "[ABC3]" "trivValue_congr(自明化の取り替えは単元倍)"
      (.inProject "ABC3" "ABC3.Found.GenEll.trivValue_congr") 3,
    .implicitStep
      ("★★★★★残っているのは**存在**である——IsTensorOf を満たす " ++
       "LocalMetric X (A ⊗ B) を mA, mB から実際に作る段。" ++
       "★そこでは「(A⊗B)|_V は自明だが A|_V は自明でない V」が問題になり、" ++
       "両方が自明になる小さい W ∋ p へ降りて貼り合わせる必要がある。" ++
       "★★値が選択に依らないことは本ファイルの normOf_indep と compat が与える") 3,
    .implicitStep
      ("★★★★★★もう 1 つ: 制限との両立(h を W ⊆ V へ落とす欄)を構造体に足す段。" ++
       "本ファイルは同じ V の上での自明化の取り替えだけを扱う") 3 ]

end ABC3.Found.Arakelov
