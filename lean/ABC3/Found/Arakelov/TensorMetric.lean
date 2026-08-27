/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.LocalMetric
import ABC3.Found.Arakelov.ArcCoverPou
import ABC3.Meta.Claim

/-!
# テンソル積の計量の**存在**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★台帳 `arakelov-coherent-base-metric` の**到達点**

| 段 | 内容 | 場所 |
|---|---|---|
| 1 | `trivValue` がテンソル積で掛け算になる | `TrivTensor.lean` |
| 2 | 切断のノルムが掛け算になる（ファイバー迂回） | `TrivSecNorm.lean` |
| 3 | 基準ノルムを担ぐ構造体・チャート独立性 | `LocalMetric.lean` |
| 4 | **テンソル計量の存在** | ★★**本ファイル** |

★2026-08-28 に `Definition 1.1` の項目全体の `.src` を下げた理由は

    `TorsorMetric.base` は対象ごとの `Classical.choice` なので
    `base_{[L·M]} ≠ base_{[L]} ⊗ base_{[M]}` となり、
    群法則が**計量のテンソル積を表さない**

であった。★★**本ファイルでその塡がりが閉じる**——局所自明な前層加群については
`mA`, `mB` から `LocalMetric X (A ⊗ B)` が**構成できて**、しかも

    `‖s ⊗ t‖ = ‖s‖ · ‖t‖`   （`normOf_tensor_metric`）

が成り立つ。

## ★★★★★★★★構成の要 —— 「`A|_V` が自明でない `V`」

`(A ⊗ B)|_V` は自明だが `A|_V` は自明でない `V` がありうる
（連結スキーム上でも `A ⊗ B ≅ 𝒪` は `A ≅ 𝒪` を含意しない）。
★そこで**局所自明性から `p` を含む `W ≤ V` で両方が自明になるものを取り**、
そこで `h_A · h_B` と置いて遷移単元で `V` の自明化 `f` へ運ぶ。

★★値がチャートに依らないことは `chartH_indep` が与える
——共通の細分へ降りて（`chartH_shrink`）、そこで `transUnit_tensorTriv` と
`compat` を使う（`chartH_triv_indep`）。

## ★★★★★★構造の恒等式は 4 つとも `rfl` であった（2026-08-28 実測）

| 恒等式 | 内容 |
|---|---|
| `trivEquiv_tensorTriv` | 自明化の値がテンソルで掛け算 |
| `trivialOfLe_tensorTriv` | テンソル自明化の制限はテンソル |
| `trivialOfLe_refl` | 自分自身への制限は恒等 |
| `trivialOfLe_trans` | 制限の推移律 |

★`set_option backward.isDefEq.respectTransparency false` が要る
（`tools\lean-idioms.md`）。

## ★残っている段（明示）

★★連続性の欄と、`APic` の群法則（同型類の群）への載せ替え。
★本ファイルは**計量そのもののテンソル積**を与えるところまでである。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite
open ABC3.Found.GenEll

variable {X : Scheme.{0}}

set_option backward.isDefEq.respectTransparency false in
theorem trivEquiv_tensorTriv {A B : X.PresheafOfModules} {V : X.Opens}
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V))
    (x : A.obj (op V)) (y : B.obj (op V)) :
    trivEquiv (A ⊗ B) V (tensorTriv eA eB) (x ⊗ₜ[(Γ(X, V) : Type)] y)
      = trivEquiv A V eA x * trivEquiv B V eB y := rfl

set_option backward.isDefEq.respectTransparency false in
theorem trivialOfLe_tensorTriv {A B : X.PresheafOfModules} {V W : X.Opens} (hWV : W ≤ V)
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V)) :
    trivialOfLe (A ⊗ B) hWV (tensorTriv eA eB)
      = tensorTriv (trivialOfLe A hWV eA) (trivialOfLe B hWV eB) := rfl

/-- ★★テンソル自明化の生成切断は生成切断のテンソルである。 -/
theorem tensorTriv_gen {A B : X.PresheafOfModules} {V : X.Opens}
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V)) :
    (trivEquiv (A ⊗ B) V (tensorTriv eA eB)).symm 1
      = ((trivEquiv A V eA).symm 1) ⊗ₜ[(Γ(X, V) : Type)] ((trivEquiv B V eB).symm 1) := by
  apply (trivEquiv (A ⊗ B) V (tensorTriv eA eB)).injective
  rw [LinearEquiv.apply_symm_apply, trivEquiv_tensorTriv,
    LinearEquiv.apply_symm_apply, LinearEquiv.apply_symm_apply, one_mul]

/-- ★★★★★★**遷移単元はテンソル積で掛け算になる**。

★これがテンソル計量の存在の**核**である——`h` の候補値が
`eA`, `eB` の取り方に依らないことは、この式と `compat` から出る。 -/
theorem transUnit_tensorTriv {A B : X.PresheafOfModules} {V : X.Opens}
    (eA eA' : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB eB' : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V)) :
    transUnit (A ⊗ B) V (tensorTriv eA eB) (tensorTriv eA' eB')
      = transUnit A V eA eA' * transUnit B V eB eB' := by
  show trivEquiv (A ⊗ B) V (tensorTriv eA' eB')
      ((trivEquiv (A ⊗ B) V (tensorTriv eA eB)).symm 1) = _
  rw [tensorTriv_gen eA eB, trivEquiv_tensorTriv]
  rfl

set_option backward.isDefEq.respectTransparency false in
/-- ★自明化を自分自身へ制限しても変わらない。 -/
theorem trivialOfLe_refl (F : X.PresheafOfModules) {V : X.Opens}
    (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V)) :
    trivialOfLe F (le_refl V) e = e := rfl

set_option backward.isDefEq.respectTransparency false in
/-- ★★制限の推移律。 -/
theorem trivialOfLe_trans (F : X.PresheafOfModules) {V W W' : X.Opens}
    (hW'W : W' ≤ W) (hWV : W ≤ V)
    (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V)) :
    trivialOfLe F (hW'W.trans hWV) e = trivialOfLe F hW'W (trivialOfLe F hWV e) := rfl

/-! ## ★★★★両方が自明になるチャート -/

/-- ★**`p` の周りで `A` と `B` の両方が自明になる `V` の中の開集合**。 -/
structure TensorChart (A B : X.PresheafOfModules) (V : X.Opens)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) where
  /-- 小さい開集合。 -/
  W : X.Opens
  /-- `V` の中にある。 -/
  hWV : W ≤ V
  /-- `p` を含む。 -/
  hpW : p ⁻¹ᵁ W = ⊤
  /-- `A` の自明化。 -/
  eA : (restrictPresheafFunctor X W).obj A ≅ 𝟙_ (PresheafModulesOn X W)
  /-- `B` の自明化。 -/
  eB : (restrictPresheafFunctor X W).obj B ≅ 𝟙_ (PresheafModulesOn X W)

/-- ★★★**局所自明なら、`p` を含む `V` の中にチャートが取れる**。

★機構は「2 つの被覆篩の交わりは被覆篩」＋「`Spec ℂ` は 1 点」。 -/
theorem nonempty_tensorChart {A B : X.PresheafOfModules}
    (hA : IsLocallyTrivial X A) (hB : IsLocallyTrivial X B) (V : X.Opens)
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) :
    Nonempty (TensorChart A B V p) := by
  obtain ⟨SA, hSA, htrivA⟩ := hA V
  obtain ⟨SB, hSB, htrivB⟩ := hB V
  have hcov := (Opens.grothendieckTopology X).intersection_covering hSA hSB
  have hxV : p.base default ∈ V := by
    have : (default : (Spec (CommRingCat.of ℂ) : Type)) ∈ p ⁻¹ᵁ V := by rw [hp]; trivial
    exact this
  obtain ⟨W, g, hg, hW⟩ := hcov (p.base default) hxV
  have hWV : W ≤ V := g.le
  have hgt : (SA ⊓ SB).arrows (homOfLE hWV) :=
    (Subsingleton.elim g (homOfLE hWV)) ▸ hg
  have htop : p ⁻¹ᵁ W = ⊤ :=
    preimage_eq_top_of_mem p W (fun z => by rw [Subsingleton.elim z default]; exact hW)
  exact ⟨⟨W, hWV, htop, (htrivA (homOfLE hWV) hgt.1).some, (htrivB (homOfLE hWV) hgt.2).some⟩⟩

/-- ★チャートを大きい開集合の上のチャートと見る。 -/
def TensorChart.lift {A B : X.PresheafOfModules} {V V' : X.Opens}
    {p : Spec (CommRingCat.of ℂ) ⟶ X} (c : TensorChart A B V p) (hVV' : V ≤ V') :
    TensorChart A B V' p :=
  ⟨c.W, c.hWV.trans hVV', c.hpW, c.eA, c.eB⟩

/-- ★チャートを小さい開集合へ降ろす。 -/
noncomputable def TensorChart.shrink {A B : X.PresheafOfModules} {V : X.Opens}
    {p : Spec (CommRingCat.of ℂ) ⟶ X} (c : TensorChart A B V p) {W' : X.Opens}
    (hW'W : W' ≤ c.W) (hpW' : p ⁻¹ᵁ W' = ⊤) : TensorChart A B V p :=
  ⟨W', hW'W.trans c.hWV, hpW', trivialOfLe A hW'W c.eA, trivialOfLe B hW'W c.eB⟩

/-! ## ★★★★★★候補値とその選択独立性 -/

/-- ★★**チャートが与える基準ノルムの候補値**。 -/
noncomputable def chartH {A B : X.PresheafOfModules} (mA : LocalMetric X A)
    (mB : LocalMetric X B) {V : X.Opens} {p : Spec (CommRingCat.of ℂ) ⟶ X}
    (c : TensorChart A B V p)
    (f : (restrictPresheafFunctor X V).obj (A ⊗ B) ≅ 𝟙_ (PresheafModulesOn X V)) : ℝ :=
  ‖evalOn p c.W c.hpW (transUnit (A ⊗ B) c.W (tensorTriv c.eA c.eB)
      (trivialOfLe (A ⊗ B) c.hWV f))‖⁻¹ * (mA.h c.W c.eA p * mB.h c.W c.eB p)

theorem chartH_pos {A B : X.PresheafOfModules} (mA : LocalMetric X A) (mB : LocalMetric X B)
    {V : X.Opens} {p : Spec (CommRingCat.of ℂ) ⟶ X} (c : TensorChart A B V p)
    (f : (restrictPresheafFunctor X V).obj (A ⊗ B) ≅ 𝟙_ (PresheafModulesOn X V)) :
    0 < chartH mA mB c f :=
  mul_pos
    (inv_pos.2 (lt_of_le_of_ne (norm_nonneg _) (Ne.symm (norm_ne_zero_iff.2
      (evalOn_ne_zero_of_isUnit p c.W c.hpW (isUnit_transUnit _ c.W _ _))))))
    (mul_pos (mA.pos c.W c.eA p) (mB.pos c.W c.eB p))

/-- ★遷移単元は互いに逆である。 -/
theorem transUnit_mul_symm (F : X.PresheafOfModules) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V)) :
    transUnit F V e e' * transUnit F V e' e = 1 :=
  (transUnit_trans F V e e' e).symm.trans (transUnit_self F V e)

/-- ★★遷移単元の値のノルムは互いに逆数である。 -/
theorem norm_evalOn_transUnit_symm (F : X.PresheafOfModules) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) :
    ‖evalOn p V hp (transUnit F V e e')‖ * ‖evalOn p V hp (transUnit F V e' e)‖ = 1 := by
  rw [← norm_mul, ← evalOn_mul, transUnit_mul_symm, evalOn_one, norm_one]

/-- ★★★★**候補値は自明化の取り方に依らない**（同じ `W` の上で）。

★機構は `transUnit_tensorTriv`（テンソルの遷移単元は積）と `compat`。 -/
theorem chartH_triv_indep {A B : X.PresheafOfModules} (mA : LocalMetric X A)
    (mB : LocalMetric X B) {V W : X.Opens} {p : Spec (CommRingCat.of ℂ) ⟶ X}
    (hWV : W ≤ V) (hpW : p ⁻¹ᵁ W = ⊤)
    (eA eA' : (restrictPresheafFunctor X W).obj A ≅ 𝟙_ (PresheafModulesOn X W))
    (eB eB' : (restrictPresheafFunctor X W).obj B ≅ 𝟙_ (PresheafModulesOn X W))
    (f : (restrictPresheafFunctor X V).obj (A ⊗ B) ≅ 𝟙_ (PresheafModulesOn X V)) :
    chartH mA mB (⟨W, hWV, hpW, eA', eB'⟩ : TensorChart A B V p) f
      = chartH mA mB (⟨W, hWV, hpW, eA, eB⟩ : TensorChart A B V p) f := by
  have hA' : mA.h W eA' p = mA.h W eA p * ‖evalOn p W hpW (transUnit A W eA' eA)‖ := by
    rw [← mA.compat W eA eA' p hpW, mul_assoc,
      norm_evalOn_transUnit_symm A W eA eA' p hpW, mul_one]
  have hB' : mB.h W eB' p = mB.h W eB p * ‖evalOn p W hpW (transUnit B W eB' eB)‖ := by
    rw [← mB.compat W eB eB' p hpW, mul_assoc,
      norm_evalOn_transUnit_symm B W eB eB' p hpW, mul_one]
  have hane : ‖evalOn p W hpW (transUnit A W eA' eA)‖ ≠ 0 :=
    norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit p W hpW (isUnit_transUnit _ W eA' eA))
  have hbne : ‖evalOn p W hpW (transUnit B W eB' eB)‖ ≠ 0 :=
    norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit p W hpW (isUnit_transUnit _ W eB' eB))
  have htne : ‖evalOn p W hpW (transUnit (A ⊗ B) W (tensorTriv eA eB)
      (trivialOfLe (A ⊗ B) hWV f))‖ ≠ 0 :=
    norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit p W hpW (isUnit_transUnit _ W _ _))
  have hfac : ‖evalOn p W hpW (transUnit (A ⊗ B) W (tensorTriv eA' eB')
        (trivialOfLe (A ⊗ B) hWV f))‖
      = ‖evalOn p W hpW (transUnit A W eA' eA)‖ * ‖evalOn p W hpW (transUnit B W eB' eB)‖
        * ‖evalOn p W hpW (transUnit (A ⊗ B) W (tensorTriv eA eB)
            (trivialOfLe (A ⊗ B) hWV f))‖ := by
    rw [transUnit_trans (A ⊗ B) W (tensorTriv eA' eB') (tensorTriv eA eB)
      (trivialOfLe (A ⊗ B) hWV f), transUnit_tensorTriv eA' eA eB' eB,
      evalOn_mul, norm_mul, evalOn_mul, norm_mul]
  show ‖evalOn p W hpW (transUnit (A ⊗ B) W (tensorTriv eA' eB')
      (trivialOfLe (A ⊗ B) hWV f))‖⁻¹ * (mA.h W eA' p * mB.h W eB' p)
    = ‖evalOn p W hpW (transUnit (A ⊗ B) W (tensorTriv eA eB)
        (trivialOfLe (A ⊗ B) hWV f))‖⁻¹ * (mA.h W eA p * mB.h W eB p)
  rw [hfac, hA', hB']
  field_simp

/-- ★★★★**候補値は小さい開集合へ降りても変わらない**。

★機構は `trivialOfLe_tensorTriv`・`trivialOfLe_trans`（ともに `rfl`）と
`transUnit_restrict`・`evalOn_restrict`・`restrict` の欄。 -/
theorem chartH_shrink {A B : X.PresheafOfModules} (mA : LocalMetric X A)
    (mB : LocalMetric X B) {V : X.Opens} {p : Spec (CommRingCat.of ℂ) ⟶ X}
    (c : TensorChart A B V p) {W' : X.Opens} (hW'W : W' ≤ c.W) (hpW' : p ⁻¹ᵁ W' = ⊤)
    (f : (restrictPresheafFunctor X V).obj (A ⊗ B) ≅ 𝟙_ (PresheafModulesOn X V)) :
    chartH mA mB (c.shrink hW'W hpW') f = chartH mA mB c f := by
  show ‖evalOn p W' hpW' (transUnit (A ⊗ B) W'
      (tensorTriv (trivialOfLe A hW'W c.eA) (trivialOfLe B hW'W c.eB))
      (trivialOfLe (A ⊗ B) (hW'W.trans c.hWV) f))‖⁻¹
      * (mA.h W' (trivialOfLe A hW'W c.eA) p * mB.h W' (trivialOfLe B hW'W c.eB) p) = _
  rw [mA.restrict hW'W c.eA p hpW', mB.restrict hW'W c.eB p hpW',
    ← trivialOfLe_tensorTriv hW'W c.eA c.eB,
    trivialOfLe_trans (A ⊗ B) hW'W c.hWV f,
    transUnit_restrict (A ⊗ B) hW'W (tensorTriv c.eA c.eB) (trivialOfLe (A ⊗ B) c.hWV f),
    evalOn_restrict p hW'W hpW']
  rfl

/-- ★★★★★★★★**候補値はチャートの取り方に全く依らない**。

★機構は「共通の細分 `c.W ⊓ c'.W` へ降りて `chartH_triv_indep`」。 -/
theorem chartH_indep {A B : X.PresheafOfModules} (mA : LocalMetric X A)
    (mB : LocalMetric X B) {V : X.Opens} {p : Spec (CommRingCat.of ℂ) ⟶ X}
    (c c' : TensorChart A B V p)
    (f : (restrictPresheafFunctor X V).obj (A ⊗ B) ≅ 𝟙_ (PresheafModulesOn X V)) :
    chartH mA mB c f = chartH mA mB c' f := by
  have hp'' : p ⁻¹ᵁ (c.W ⊓ c'.W) = ⊤ := by
    show p ⁻¹ᵁ c.W ⊓ p ⁻¹ᵁ c'.W = ⊤
    rw [c.hpW, c'.hpW, inf_idem]
  rw [← chartH_shrink mA mB c (inf_le_left : c.W ⊓ c'.W ≤ c.W) hp'' f,
    ← chartH_shrink mA mB c' (inf_le_right : c.W ⊓ c'.W ≤ c'.W) hp'' f]
  exact chartH_triv_indep mA mB _ hp'' _ _ _ _ f

/-! ## ★★★★★★★★★テンソル積の計量 -/

open scoped Classical in
/-- ★★★★★★★★★**テンソル積の計量**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★`p` の周りで両方が自明になるチャートを取り、そこで `h_A · h_B` と置いて
遷移単元で `V` の自明化 `f` へ運ぶ。★★値がチャートに依らないことは `chartH_indep`。 -/
noncomputable def LocalMetric.tensor {A B : X.PresheafOfModules}
    (hA : IsLocallyTrivial X A) (hB : IsLocallyTrivial X B)
    (mA : LocalMetric X A) (mB : LocalMetric X B) : LocalMetric X (A ⊗ B) where
  h V f p := if hc : Nonempty (TensorChart A B V p) then chartH mA mB hc.some f else 1
  pos V f p := by
    by_cases hc : Nonempty (TensorChart A B V p)
    · show 0 < (if hc : Nonempty (TensorChart A B V p) then chartH mA mB hc.some f else 1)
      rw [dif_pos hc]
      exact chartH_pos mA mB hc.some f
    · show 0 < (if hc : Nonempty (TensorChart A B V p) then chartH mA mB hc.some f else 1)
      rw [dif_neg hc]
      exact one_pos
  compat V e e' p hp := by
    have hc : Nonempty (TensorChart A B V p) := nonempty_tensorChart hA hB V p hp
    show (if hc : Nonempty (TensorChart A B V p) then chartH mA mB hc.some e' else 1) * _
      = (if hc : Nonempty (TensorChart A B V p) then chartH mA mB hc.some e else 1)
    rw [dif_pos hc, dif_pos hc]
    set c := hc.some with hcdef
    have hstep : ‖evalOn p c.W c.hpW (transUnit (A ⊗ B) c.W (tensorTriv c.eA c.eB)
          (trivialOfLe (A ⊗ B) c.hWV e'))‖
        = ‖evalOn p c.W c.hpW (transUnit (A ⊗ B) c.W (tensorTriv c.eA c.eB)
            (trivialOfLe (A ⊗ B) c.hWV e))‖ * ‖evalOn p V hp (transUnit (A ⊗ B) V e e')‖ := by
      rw [transUnit_trans (A ⊗ B) c.W (tensorTriv c.eA c.eB)
        (trivialOfLe (A ⊗ B) c.hWV e) (trivialOfLe (A ⊗ B) c.hWV e'),
        transUnit_restrict (A ⊗ B) c.hWV e e', evalOn_mul, norm_mul,
        evalOn_restrict p c.hWV c.hpW]
    have hne : ‖evalOn p c.W c.hpW (transUnit (A ⊗ B) c.W (tensorTriv c.eA c.eB)
        (trivialOfLe (A ⊗ B) c.hWV e))‖ ≠ 0 :=
      norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit p c.W c.hpW (isUnit_transUnit _ c.W _ _))
    have hune : ‖evalOn p V hp (transUnit (A ⊗ B) V e e')‖ ≠ 0 :=
      norm_ne_zero_iff.2 (evalOn_ne_zero_of_isUnit p V hp (isUnit_transUnit _ V e e'))
    show ‖evalOn p c.W c.hpW (transUnit (A ⊗ B) c.W (tensorTriv c.eA c.eB)
          (trivialOfLe (A ⊗ B) c.hWV e'))‖⁻¹ * (mA.h c.W c.eA p * mB.h c.W c.eB p)
        * ‖evalOn p V hp (transUnit (A ⊗ B) V e e')‖
      = ‖evalOn p c.W c.hpW (transUnit (A ⊗ B) c.W (tensorTriv c.eA c.eB)
          (trivialOfLe (A ⊗ B) c.hWV e))‖⁻¹ * (mA.h c.W c.eA p * mB.h c.W c.eB p)
    rw [hstep, mul_inv]
    field_simp
  restrict {V W} hWV e p hpW := by
    have hcW : Nonempty (TensorChart A B W p) := nonempty_tensorChart hA hB W p hpW
    have hcV : Nonempty (TensorChart A B V p) :=
      nonempty_tensorChart hA hB V p (preimage_eq_top_of_le hWV hpW)
    show (if hc : Nonempty (TensorChart A B W p) then
        chartH mA mB hc.some (trivialOfLe (A ⊗ B) hWV e) else 1)
      = (if hc : Nonempty (TensorChart A B V p) then chartH mA mB hc.some e else 1)
    rw [dif_pos hcW, dif_pos hcV]
    have hlift : chartH mA mB (hcW.some.lift hWV) e
        = chartH mA mB hcW.some (trivialOfLe (A ⊗ B) hWV e) := rfl
    rw [← hlift]
    exact chartH_indep mA mB _ _ e

/-- ★★★★★★★★★**それは実際にテンソル積の計量である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★これで `normOf_tensor`（`‖s ⊗ t‖ = ‖s‖ · ‖t‖`）が**無条件に**使える。 -/
theorem isTensorOf_tensor {A B : X.PresheafOfModules}
    (hA : IsLocallyTrivial X A) (hB : IsLocallyTrivial X B)
    (mA : LocalMetric X A) (mB : LocalMetric X B) :
    IsTensorOf mA mB (LocalMetric.tensor hA hB mA mB) := by
  classical
  intro V eA eB p hp
  have hc : Nonempty (TensorChart A B V p) := nonempty_tensorChart hA hB V p hp
  show (if hc : Nonempty (TensorChart A B V p) then
      chartH mA mB hc.some (tensorTriv eA eB) else 1) = _
  rw [dif_pos hc]
  have hself : chartH mA mB (⟨V, le_refl V, hp, eA, eB⟩ : TensorChart A B V p)
      (tensorTriv eA eB) = mA.h V eA p * mB.h V eB p := by
    show ‖evalOn p V hp (transUnit (A ⊗ B) V (tensorTriv eA eB)
        (trivialOfLe (A ⊗ B) (le_refl V) (tensorTriv eA eB)))‖⁻¹ * _ = _
    rw [trivialOfLe_refl (A ⊗ B) (tensorTriv eA eB), transUnit_self, evalOn_one,
      norm_one, inv_one, one_mul]
  rw [← hself]
  exact chartH_indep mA mB _ _ (tensorTriv eA eB)

/-- ★★★★★★★★★**`‖s ⊗ t‖ = ‖s‖ · ‖t‖`——無条件の形**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★これが台帳 `arakelov-coherent-base-metric` の**到達点**である
——`Classical.choice` の基準計量では落ちていた加法性が、
局所自明な前層加群の計量については**構成つきで**出る。 -/
theorem normOf_tensor_metric {A B : X.PresheafOfModules}
    (hA : IsLocallyTrivial X A) (hB : IsLocallyTrivial X B)
    (mA : LocalMetric X A) (mB : LocalMetric X B) (V : X.Opens)
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤)
    (s : A.obj (op ⊤)) (t : B.obj (op ⊤)) :
    (LocalMetric.tensor hA hB mA mB).normOf V (tensorTriv eA eB) p hp
        (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t)
      = mA.normOf V eA p hp s * mB.normOf V eB p hp t :=
  normOf_tensor (isTensorOf_tensor hA hB mA mB) V eA eB p hp s t

/-- ★★★★★★★**テンソル自明化を選ばない形**——任意の自明化 `f` で見ても掛け算になる。 -/
theorem normOf_tensor_metric' {A B : X.PresheafOfModules}
    (hA : IsLocallyTrivial X A) (hB : IsLocallyTrivial X B)
    (mA : LocalMetric X A) (mB : LocalMetric X B) (V : X.Opens)
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V))
    (f : (restrictPresheafFunctor X V).obj (A ⊗ B) ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤)
    (s : A.obj (op ⊤)) (t : B.obj (op ⊤)) :
    (LocalMetric.tensor hA hB mA mB).normOf V f p hp
        (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t)
      = mA.normOf V eA p hp s * mB.normOf V eB p hp t :=
  normOf_tensor' (isTensorOf_tensor hA hB mA mB) V eA eB f p hp s t

/-! ### ★出典の紐付け(`.src`) -/

def transUnit_tensorTriv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(遷移単元はテンソル積で掛け算になること)",
    sectionId := "genell-def-1-1-i" }

def nonempty_tensorChart.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(局所自明なら両方が自明になるチャートが取れること)",
    sectionId := "genell-def-1-1-i" }

def chartH_indep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(テンソル計量の候補値はチャートに依らないこと)",
    sectionId := "genell-def-1-1-i" }

def LocalMetric.tensor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(テンソル積の計量の構成)",
    sectionId := "genell-def-1-1-i" }

def isTensorOf_tensor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(構成した計量が実際にテンソル積の計量であること)",
    sectionId := "genell-def-1-1-i" }

def normOf_tensor_metric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(‖s ⊗ t‖ = ‖s‖ · ‖t‖——無条件の形)",
    sectionId := "genell-def-1-1-i" }

def LocalMetric.tensor.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "trivValue_tensor(自明化の値がテンソル積で掛け算になる)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.trivValue_tensor") 3,
    .citation "[ABC3]" "transUnit_restrict / evalOn_restrict(制限との両立)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.transUnit_restrict") 3,
    .citation "[ABC3]" "IsLocallyTrivial(可逆層の強い局所自明性)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.IsLocallyTrivial") 3,
    .implicitStep
      ("★★★★★構成の要は「(A⊗B)|_V は自明だが A|_V は自明でない V」の扱いである。" ++
       "★局所自明性から p を含む W ≤ V で両方が自明になるものを取り、" ++
       "そこで h_A · h_B と置いて遷移単元で V の自明化 f へ運ぶ。" ++
       "★★値がチャートに依らないことは chartH_indep(共通の細分へ降りて " ++
       "transUnit_tensorTriv と compat)が与える") 3,
    .implicitStep
      ("★★★★★★構造の恒等式は 4 つとも rfl であった(2026-08-28 実測): " ++
       "trivValue_tensor・trivEquiv_tensorTriv・trivialOfLe_tensorTriv・" ++
       "trivialOfLe_refl / trivialOfLe_trans。" ++
       "set_option backward.isDefEq.respectTransparency false が要る") 3,
    .implicitStep
      ("★★残る段: 連続性の欄と、APic の群法則(同型類の群)への載せ替え。" ++
       "★本ファイルは計量そのもののテンソル積を与えるところまでである") 3 ]

end ABC3.Found.Arakelov
