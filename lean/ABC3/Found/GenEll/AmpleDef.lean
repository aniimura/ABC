/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.PicLocalTrivial
import ABC3.Found.Arakelov.PicBasicTrivial
import Mathlib.AlgebraicGeometry.AffineScheme
import ABC3.Meta.Claim

/-!
# 切断の**非消失軌跡 `X_s`** と **ample の定義**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★`ample-and-projective-embedding` の段 D

`ample` の定義（[Stacks] 28.27.1）は

> `X` が準コンパクトで、各点 `x` に `n ≥ 1` と `s ∈ Γ(X, L^{⊗n})` があって
> `x ∈ X_s` かつ `X_s` がアフィンであること

である。★核は **`X_s`（切断の非消失軌跡）** であり、本ファイルはそれを作る。

## ★★★★★★mathlib に加群層の**茎**は無い —— 自明化で回る

`X_s` の教科書的な定義は「茎 `M_x` の中で `s_x ∉ 𝔪_x M_x`」だが、
★**mathlib に `SheafOfModules` の茎は 1 件も無い**（2026-08-27 実測）。

★★そこで**局所自明化を使う**。`IsLocallyTrivial`（`Found/Arakelov/PicLocalTrivial.lean`）が
与える `e : M|_V ≅ 𝟙_` で切断を関数 `f_V ∈ Γ(V, 𝒪_X)` と見て、
`Scheme.basicOpen f_V` を取り、**すべての `(V, e)` について合併する**:

    `X_s ≝ ⨆_{V} ⨆_{e : M|_V ≅ 𝟙_} X.basicOpen (trivValue M V e s)`

★★★これは**選択を使わない**（自明化を選ばない）ので、
定義そのものは `IsLocallyTrivial` を仮定しなくても書ける。開であることは合併だから自動。

## ★★★★★自明化の取り方に依らないこと

同じ `V` の上の 2 つの自明化 `e, e'` については

    `trivValue e' s = trivValue e s · u`（`u` は単元）

が成り立つ（`trivValue_congr`）。★機構は「環 `R` を自分自身の上の加群と見たときの
線型同型は `θ 1` 倍であり、`θ 1` は単元」という**純代数**である
（`linearEquiv_self_apply` / `linearEquiv_self_isUnit`）。
★★したがって `basicOpen` は変わらない（`basicOpen_trivValue_congr`）。

## ★★★★★★★異なる `V, W` にまたがる独立性（段 D2）

    `X_s ⊓ V = X.basicOpen (trivValue M V e s)`（`nonVanishing_inf`）

★機構は `trivialOfLe`（`Found/Arakelov/PicBasicTrivial.lean`）と `trivValue` が
**可換**であること（`trivValue_restrict`）だけで、それは `e.hom` の**自然性**から出る。
★★`trivialOfLe` は 2 つの `rfl` 等式を `rw` して作られているので、
**評価すると計算する**——証明中の `h1`・`h2` はどちらも `rfl` である。

## ★★★★★★定義の健全性

構造層 `𝒪_X` については **`X_s` は `Scheme.basicOpen s` そのもの**である
（`nonVanishing_unit`）。★これで `nonVanishing` が
「知っている `basicOpen` の一般化」だと型で言える——**空虚ではない**。

## ★残っている段（明示）

★**very ample（`ℙᴺ` への閉埋め込み）**は段 C2（`O(1)`）の後である。
★★**Serre の定理**（ample のある冪が very ample）は段 E。
★★★**座標の高さと `htArith` の比較**がその先にある。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

universe u

/-! ## ★環を自分自身の上の加群と見たときの線型同型 -/

/-- ★**`θ f = f · θ 1`** —— `R` を `R`-加群と見たときの線型写像は `θ 1` 倍である。 -/
theorem linearEquiv_self_apply {R : Type*} [Ring R] (θ : R ≃ₗ[R] R) (f : R) :
    θ f = f * θ 1 := by
  have : θ (f • (1 : R)) = f • θ 1 := θ.map_smul f 1
  simpa [smul_eq_mul] using this

/-- ★★**`θ 1` は単元**——逆写像も `θ⁻¹ 1` 倍だからである。 -/
theorem linearEquiv_self_isUnit {R : Type*} [Ring R] (θ : R ≃ₗ[R] R) : IsUnit (θ 1) := by
  have h1 : θ.symm (θ 1) = 1 := θ.symm_apply_apply 1
  rw [linearEquiv_self_apply θ.symm (θ 1)] at h1
  have h2 : θ (θ.symm 1) = 1 := θ.apply_symm_apply 1
  rw [linearEquiv_self_apply θ (θ.symm 1)] at h2
  exact ⟨⟨θ 1, θ.symm 1, h1, h2⟩, rfl⟩

variable {X : Scheme.{u}}

/-! ## ★★★切断を自明化で関数と見る -/

/-- ★大域切断を開集合 `V` へ制限する。 -/
noncomputable def secOn (M : X.PresheafOfModules) (V : X.Opens)
    (s : M.obj (op ⊤)) : M.obj (op V) :=
  M.map (homOfLE le_top : V ⟶ (⊤ : X.Opens)).op s

/-- ★★**自明化 `e : M|_V ≅ 𝟙_` を `Γ(V, 𝒪_X)` 値の線型同型として読む**。

★`Over V` の終対象 `⟨V, 𝟙⟩` で評価するだけである
（`IsLocallyTrivial.isLocallyRankOne` と同じ機構）。 -/
noncomputable def trivEquiv (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V)) :
    (M.obj (op V) : Type u) ≃ₗ[(Γ(X, V) : Type u)] (Γ(X, V) : Type u) :=
  ((PresheafOfModules.evaluation _ (op (Over.mk (𝟙 V)))).mapIso e).toLinearEquiv

/-- ★★★**自明化 `e` のもとで切断 `s` が定める関数** `f_V ∈ Γ(V, 𝒪_X)`。 -/
noncomputable def trivValue (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s : M.obj (op ⊤)) : Γ(X, V) :=
  ((PresheafOfModules.evaluation _ (op (Over.mk (𝟙 V)))).mapIso e).hom (secOn M V s)

theorem trivValue_eq_trivEquiv (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s : M.obj (op ⊤)) :
    trivValue M V e s = trivEquiv M V e (secOn M V s) := rfl

theorem trivValue_zero (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V)) :
    trivValue M V e (0 : M.obj (op ⊤)) = 0 := by
  simp only [trivValue, secOn, map_zero]
  exact map_zero _

/-! ## ★★★★★自明化の取り方に依らないこと（同じ `V` の上で） -/

/-- ★★★★★**同じ `V` の上の 2 つの自明化は単元倍で違う**。

★機構は `linearEquiv_self_apply` / `linearEquiv_self_isUnit`（純代数）だけである。 -/
theorem trivValue_congr (M : X.PresheafOfModules) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s : M.obj (op ⊤)) :
    ∃ u : Γ(X, V), IsUnit u ∧ trivValue M V e' s = trivValue M V e s * u := by
  set θ : (Γ(X, V) : Type u) ≃ₗ[(Γ(X, V) : Type u)] (Γ(X, V) : Type u) :=
    (trivEquiv M V e).symm.trans (trivEquiv M V e') with hθ
  refine ⟨θ 1, linearEquiv_self_isUnit θ, ?_⟩
  have h := linearEquiv_self_apply θ (trivEquiv M V e (secOn M V s))
  have h2 : θ (trivEquiv M V e (secOn M V s)) = trivEquiv M V e' (secOn M V s) := by
    rw [hθ]; simp
  rw [trivValue_eq_trivEquiv, trivValue_eq_trivEquiv, ← h2, h]

/-- ★★★★★★**したがって `basicOpen` は自明化の取り方に依らない**（同じ `V` の上で）。 -/
theorem basicOpen_trivValue_congr (M : X.PresheafOfModules) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s : M.obj (op ⊤)) :
    X.basicOpen (trivValue M V e' s) = X.basicOpen (trivValue M V e s) := by
  obtain ⟨u, hu, huv⟩ := trivValue_congr M V e e' s
  rw [huv, X.basicOpen_mul, X.basicOpen_of_isUnit hu, inf_eq_left]
  exact X.basicOpen_le _

/-- ★★★★★**単元は切断に依らない**（`trivValue_congr` の強い形）。

★★★これが「切断の比 `s/t` が自明化に依らない」ことの根拠になる
——分子と分母が**同じ** `u` 倍されるからである。 -/
theorem trivValue_congr' (M : X.PresheafOfModules) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V)) :
    ∃ u : Γ(X, V), IsUnit u ∧
      ∀ s : M.obj (op ⊤), trivValue M V e' s = trivValue M V e s * u := by
  set θ : (Γ(X, V) : Type u) ≃ₗ[(Γ(X, V) : Type u)] (Γ(X, V) : Type u) :=
    (trivEquiv M V e).symm.trans (trivEquiv M V e') with hθ
  refine ⟨θ 1, linearEquiv_self_isUnit θ, fun s => ?_⟩
  have h := linearEquiv_self_apply θ (trivEquiv M V e (secOn M V s))
  have h2 : θ (trivEquiv M V e (secOn M V s)) = trivEquiv M V e' (secOn M V s) := by
    rw [hθ]; simp
  rw [trivValue_eq_trivEquiv, trivValue_eq_trivEquiv, ← h2, h]

/-! ## ★★★★★★★★非消失軌跡 `X_s` -/

/-- ★★★★★★★★**切断 `s` の非消失軌跡 `X_s`**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★**開集合として定義する**——合併だから開であることは自動である。
★★自明化を**選ばない**ので選択公理を使わない。 -/
noncomputable def nonVanishing (M : X.PresheafOfModules) (s : M.obj (op ⊤)) : X.Opens :=
  ⨆ (V : X.Opens), ⨆ (e : (restrictPresheafFunctor X V).obj M
      ≅ 𝟙_ (PresheafModulesOn X V)), X.basicOpen (trivValue M V e s)

theorem mem_nonVanishing_iff (M : X.PresheafOfModules) (s : M.obj (op ⊤)) (x : X) :
    x ∈ nonVanishing M s ↔ ∃ (V : X.Opens) (e : (restrictPresheafFunctor X V).obj M
      ≅ 𝟙_ (PresheafModulesOn X V)), x ∈ X.basicOpen (trivValue M V e s) := by
  simp [nonVanishing]

theorem basicOpen_trivValue_le (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s : M.obj (op ⊤)) :
    X.basicOpen (trivValue M V e s) ≤ nonVanishing M s :=
  le_iSup_of_le V (le_iSup_of_le e le_rfl)

/-- ★★★**零切断はどこでも消える**——空虚封じ。 -/
theorem nonVanishing_zero (M : X.PresheafOfModules) :
    nonVanishing M (0 : M.obj (op ⊤)) = ⊥ := by
  simp [nonVanishing, trivValue_zero]

/-! ## ★★★★テンソル冪 -/

/-- ★★**前層加群のテンソル冪 `M^{⊗n}`**（右結合の再帰）。

★★★mathlib の `PresheafOfModules.monoidalCategory` は**結合子・単位子まで完成している**ので、
層の水準（`Found/GenEll/InvertibleSheaf.lean`）と違って回避経路ではない。 -/
noncomputable def presheafTensorPow (M : X.PresheafOfModules) : ℕ → X.PresheafOfModules
  | 0 => 𝟙_ X.PresheafOfModules
  | n + 1 => M ⊗ presheafTensorPow M n

@[simp] theorem presheafTensorPow_zero (M : X.PresheafOfModules) :
    presheafTensorPow M 0 = 𝟙_ X.PresheafOfModules := rfl

@[simp] theorem presheafTensorPow_succ (M : X.PresheafOfModules) (n : ℕ) :
    presheafTensorPow M (n + 1) = M ⊗ presheafTensorPow M n := rfl

/-- ★構造層は局所自明である。 -/
theorem isLocallyTrivial_unit : IsLocallyTrivial X (𝟙_ X.PresheafOfModules) := by
  intro U
  refine ⟨⊤, (Opens.grothendieckTopology X).top_mem U, ?_⟩
  intro V i _
  exact ⟨restrictPresheafUnit.symm⟩

/-- ★★**局所自明性はテンソル冪で保たれる**——`IsLocallyTrivial.tensor` の帰納。 -/
theorem isLocallyTrivial_presheafTensorPow {M : X.PresheafOfModules}
    (hM : IsLocallyTrivial X M) : ∀ n : ℕ, IsLocallyTrivial X (presheafTensorPow M n)
  | 0 => isLocallyTrivial_unit
  | n + 1 => hM.tensor (isLocallyTrivial_presheafTensorPow hM n)

/-! ## ★★★★★★★★★異なる `V, W` にまたがる独立性（段 D2） -/

/-- ★**大域切断の制限は推移的**——`(s|_V)|_W = s|_W`。 -/
theorem secOn_restrict (M : X.PresheafOfModules) {V W : X.Opens} (h : W ≤ V)
    (s : M.obj (op ⊤)) :
    ((restrictPresheafFunctor X V).obj M).map
        (Over.homMk (homOfLE h) : Over.mk (homOfLE h : W ⟶ V) ⟶ Over.mk (𝟙 V)).op
        (secOn M V s) = secOn M W s := by
  simp only [secOn]
  show (ConcreteCategory.hom (M.map (homOfLE h).op)) _ = _
  rw [← PresheafOfModules.map_comp_apply]
  rfl

/-- ★★★★★★★**自明化の制限と `trivValue` は可換**。

    `trivValue M W (trivialOfLe M h e) s = (trivValue M V e s)|_W`

★機構は `e.hom` の**自然性**（`PresheafOfModules.naturality_apply`）だけである。
★★`trivialOfLe`（`Found/Arakelov/PicBasicTrivial.lean`）は 2 つの `rfl` 等式を
`rw` して作られているので、**評価すると計算する**——`h1`・`h2` はどちらも `rfl` である。 -/
theorem trivValue_restrict (M : X.PresheafOfModules) {V W : X.Opens} (h : W ≤ V)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s : M.obj (op ⊤)) :
    trivValue M W (trivialOfLe M h e) s
      = X.presheaf.map (homOfLE h).op (trivValue M V e s) := by
  have hnat := PresheafOfModules.naturality_apply e.hom
    (Y := op (Over.mk (homOfLE h : W ⟶ V)))
    (X := op (Over.mk (𝟙 V)))
    (Over.homMk (homOfLE h) : Over.mk (homOfLE h : W ⟶ V) ⟶ Over.mk (𝟙 V)).op
    (secOn M V s)
  rw [secOn_restrict M h s] at hnat
  have h1 : trivValue M W (trivialOfLe M h e) s
      = e.hom.app (op (Over.mk (homOfLE h : W ⟶ V))) (secOn M W s) := rfl
  have h2 : X.presheaf.map (homOfLE h).op (trivValue M V e s)
      = (𝟙_ (PresheafModulesOn X V)).map
          (Over.homMk (homOfLE h) : Over.mk (homOfLE h : W ⟶ V) ⟶ Over.mk (𝟙 V)).op
          (e.hom.app (op (Over.mk (𝟙 V))) (secOn M V s)) := rfl
  rw [h1, h2]
  exact hnat

theorem basicOpen_trivValue_restrict (M : X.PresheafOfModules) {V W : X.Opens} (h : W ≤ V)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s : M.obj (op ⊤)) :
    X.basicOpen (trivValue M W (trivialOfLe M h e) s)
      = W ⊓ X.basicOpen (trivValue M V e s) := by
  rw [trivValue_restrict M h e s, X.basicOpen_res]

/-- ★★★★★★★★★**`X_s` は自明化の上でちょうど `basicOpen` である**（段 D2 の到達点）。

    `X_s ⊓ V = X.basicOpen (trivValue M V e s)`

★これで `X_s` が「どの `(V, e)` から見ても同じもの」であることが型で言える。
★★機構は `V ⊓ W` へ両方の自明化を落として `basicOpen_trivValue_congr` を当てるだけ。 -/
theorem nonVanishing_inf (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s : M.obj (op ⊤)) :
    nonVanishing M s ⊓ V = X.basicOpen (trivValue M V e s) := by
  refine le_antisymm ?_ (le_inf (basicOpen_trivValue_le M V e s) (X.basicOpen_le _))
  rintro x ⟨hxN, hxV⟩
  obtain ⟨W, e', hxW⟩ := (mem_nonVanishing_iff M s x).1 hxN
  have hxW' : x ∈ W := X.basicOpen_le (trivValue M W e' s) hxW
  have hxU : x ∈ V ⊓ W := ⟨hxV, hxW'⟩
  have h1 : x ∈ X.basicOpen (trivValue M (V ⊓ W) (trivialOfLe M inf_le_right e') s) := by
    rw [basicOpen_trivValue_restrict M (inf_le_right : V ⊓ W ≤ W) e' s]
    exact ⟨hxU, hxW⟩
  rw [basicOpen_trivValue_congr M (V ⊓ W) (trivialOfLe M inf_le_left e)
    (trivialOfLe M inf_le_right e') s] at h1
  rw [basicOpen_trivValue_restrict M (inf_le_left : V ⊓ W ≤ V) e s] at h1
  exact h1.2

/-- ★★★★**局所判定**——自明化のある近傍の点については `basicOpen` を見ればよい。 -/
theorem mem_nonVanishing_iff_of_mem (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s : M.obj (op ⊤)) {x : X} (hx : x ∈ V) :
    x ∈ nonVanishing M s ↔ x ∈ X.basicOpen (trivValue M V e s) := by
  constructor
  · intro hxN
    have hxi : x ∈ nonVanishing M s ⊓ V := ⟨hxN, hx⟩
    rwa [nonVanishing_inf M V e s] at hxi
  · exact fun hxB => basicOpen_trivValue_le M V e s hxB

set_option backward.isDefEq.respectTransparency false in
theorem trivValue_unit_top (s : (𝟙_ X.PresheafOfModules).obj (op (⊤ : X.Opens))) :
    trivValue (𝟙_ X.PresheafOfModules) ⊤ restrictPresheafUnit.symm s = s := by
  simp [trivValue, secOn, restrictPresheafUnit, Functor.Monoidal.εIso]
  rfl

/-- ★★★★★★★★**定義の健全性** —— 構造層 `𝒪_X` については
`X_s` は `Scheme.basicOpen s` そのものである。

★★これで `nonVanishing` が「知っている `basicOpen` の一般化」だと型で言える
——空虚ではない。 -/
theorem nonVanishing_unit (s : (𝟙_ X.PresheafOfModules).obj (op (⊤ : X.Opens))) :
    nonVanishing (𝟙_ X.PresheafOfModules) s = X.basicOpen s := by
  have h := nonVanishing_inf (𝟙_ X.PresheafOfModules) ⊤ restrictPresheafUnit.symm s
  rw [inf_top_eq, trivValue_unit_top] at h
  exact h

/-! ## ★★★★★★★`ample` の定義 -/

/-- ★★★★★★★**ample な可逆層**（[Stacks] 28.27.1）。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

`X` が準コンパクトで、各点 `x` に `n ≥ 1` と `s ∈ Γ(X, M^{⊗n})` があって
`x ∈ X_s` かつ `X_s` がアフィンであること。

★★**可逆性（`IsLocallyTrivial`）は定義に含めない**——[Stacks] もそうしている
（ample は `𝒪_X`-加群層一般に対する述語で、そこから可逆性が従うわけではない）。
消費側では `IsLocallyTrivial` と併せて使う。 -/
def IsAmple (M : X.PresheafOfModules) : Prop :=
  CompactSpace X ∧
    ∀ x : X, ∃ n : ℕ, 0 < n ∧ ∃ s : (presheafTensorPow M n).obj (op ⊤),
      x ∈ nonVanishing (presheafTensorPow M n) s ∧
        IsAffineOpen (nonVanishing (presheafTensorPow M n) s)

/-! ## ★★★★★★★★切断の比 `s/t` —— 同次座標の入口 -/

/-- ★`basicOpen f` の中の開集合の上では `f` の制限は単元である。 -/
theorem isUnit_res_of_le_basicOpen {V U : X.Opens} (f : Γ(X, V)) (h : U ≤ V)
    (hU : U ≤ X.basicOpen f) :
    IsUnit (X.presheaf.map (homOfLE h).op f) := by
  have h1 := X.toRingedSpace.isUnit_res_basicOpen f
  have h2 := h1.map (X.presheaf.map (homOfLE hU).op).hom
  have h3 : (X.presheaf.map (homOfLE hU).op)
      ((X.presheaf.map (homOfLE (X.basicOpen_le f)).op) f)
      = X.presheaf.map (homOfLE h).op f := by
    rw [← CommRingCat.comp_apply, ← Functor.map_comp]
    rfl
  rw [← h3]
  exact h2

/-- ★★**分母は `X_t ⊓ V` の上で単元**——`nonVanishing_inf` がそのまま効く。 -/
theorem isUnit_trivValue_res (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (t : M.obj (op ⊤)) :
    IsUnit (X.presheaf.map (homOfLE (inf_le_right : nonVanishing M t ⊓ V ≤ V)).op
      (trivValue M V e t)) :=
  isUnit_res_of_le_basicOpen (trivValue M V e t) inf_le_right
    (le_of_eq (nonVanishing_inf M V e t))

/-- ★★★★★★★★**切断の比 `s/t`**——`X_t ⊓ V` の上の正則関数。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★**型が `e` を含まない**のが要点である——`nonVanishing_inf` によって
`X.basicOpen (trivValue M V e t)` を `X_t ⊓ V` と書き換えられるからで、
これで「自明化に依らない」が**転送なしの等式**として言える（`sectionRatio_congr`）。

★★★これが `northcott_of_projModel` が要求する**同次座標**の入口である
——原文の『射影埋め込み』の座標は `s_i / s_j` に他ならない。 -/
noncomputable def sectionRatio (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s t : M.obj (op ⊤)) : Γ(X, nonVanishing M t ⊓ V) :=
  X.presheaf.map (homOfLE (inf_le_right : nonVanishing M t ⊓ V ≤ V)).op (trivValue M V e s)
    * ↑(isUnit_trivValue_res M V e t).unit⁻¹

/-- ★★★**比の特徴づけ**: `(s/t) · t = s`。 -/
theorem sectionRatio_mul (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s t : M.obj (op ⊤)) :
    sectionRatio M V e s t
        * X.presheaf.map (homOfLE (inf_le_right : nonVanishing M t ⊓ V ≤ V)).op
            (trivValue M V e t)
      = X.presheaf.map (homOfLE (inf_le_right : nonVanishing M t ⊓ V ≤ V)).op
          (trivValue M V e s) := by
  have hu := isUnit_trivValue_res M V e t
  have key : (↑(hu.unit⁻¹) : Γ(X, nonVanishing M t ⊓ V)) * (↑hu.unit : Γ(X, nonVanishing M t ⊓ V))
      = 1 := hu.unit.inv_mul
  rw [hu.unit_spec] at key
  rw [sectionRatio, mul_assoc, key, mul_one]

/-- ★★★★**その特徴づけは一意である**（分母が単元だから）。 -/
theorem sectionRatio_unique (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s t : M.obj (op ⊤)) (r : Γ(X, nonVanishing M t ⊓ V))
    (hr : r * X.presheaf.map (homOfLE (inf_le_right : nonVanishing M t ⊓ V ≤ V)).op
            (trivValue M V e t)
        = X.presheaf.map (homOfLE (inf_le_right : nonVanishing M t ⊓ V ≤ V)).op
            (trivValue M V e s)) :
    r = sectionRatio M V e s t := by
  have hu := isUnit_trivValue_res M V e t
  have key : (↑hu.unit : Γ(X, nonVanishing M t ⊓ V)) * (↑(hu.unit⁻¹) : Γ(X, nonVanishing M t ⊓ V))
      = 1 := hu.unit.mul_inv
  rw [hu.unit_spec] at key
  calc r = r * (X.presheaf.map (homOfLE (inf_le_right : nonVanishing M t ⊓ V ≤ V)).op
              (trivValue M V e t) * ↑(hu.unit⁻¹)) := by rw [key, mul_one]
    _ = (r * X.presheaf.map (homOfLE (inf_le_right : nonVanishing M t ⊓ V ≤ V)).op
              (trivValue M V e t)) * ↑(hu.unit⁻¹) := by rw [mul_assoc]
    _ = sectionRatio M V e s t := by rw [hr]; rfl

/-- ★★★★★★★★★**切断の比は自明化の取り方に依らない**。

★機構は `trivValue_congr'`——**分子と分母が同じ単元 `u` 倍される**ので比では消える。
★★これで `s_i / s_j` が「`X_{s_j}` の上の正則関数」として意味を持つ。 -/
theorem sectionRatio_congr (M : X.PresheafOfModules) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (s t : M.obj (op ⊤)) :
    sectionRatio M V e' s t = sectionRatio M V e s t := by
  obtain ⟨u, hu, hall⟩ := trivValue_congr' M V e e'
  set ρ := X.presheaf.map (homOfLE (inf_le_right : nonVanishing M t ⊓ V ≤ V)).op with hρ
  have huρ : IsUnit (ρ u) := hu.map (X.presheaf.map
    (homOfLE (inf_le_right : nonVanishing M t ⊓ V ≤ V)).op).hom
  refine sectionRatio_unique M V e s t (sectionRatio M V e' s t) ?_
  have h1 := sectionRatio_mul M V e' s t
  rw [hall t, hall s, map_mul, map_mul] at h1
  have h2 : sectionRatio M V e' s t * ρ (trivValue M V e t) * ρ u
      = ρ (trivValue M V e s) * ρ u := by
    rw [mul_assoc]; exact h1
  exact IsUnit.mul_right_cancel huρ h2

@[simp] theorem sectionRatio_self (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (t : M.obj (op ⊤)) :
    sectionRatio M V e t t = 1 :=
  (sectionRatio_unique M V e t t 1 (by rw [one_mul])).symm

@[simp] theorem sectionRatio_zero (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (t : M.obj (op ⊤)) :
    sectionRatio M V e 0 t = 0 :=
  (sectionRatio_unique M V e 0 t 0 (by rw [zero_mul, trivValue_zero, map_zero])).symm

/-! ## ★★★★★★★同型に沿った移送と、`IsAmple` の空虚封じ -/

/-- ★大域切断は射に沿って移る（自然性）。 -/
theorem secOn_map (M N : X.PresheafOfModules) (φ : M ⟶ N) (V : X.Opens) (s : M.obj (op ⊤)) :
    secOn N V (φ.app (op ⊤) s) = φ.app (op V) (secOn M V s) :=
  (PresheafOfModules.naturality_apply φ (homOfLE (le_top : V ≤ ⊤)).op s).symm

/-- ★★**同型 `i : M ≅ N` は自明化を移す**——
`trivValue N V e (i s) = trivValue M V (i|_V ≪≫ e) s`。 -/
theorem trivValue_of_iso (M N : X.PresheafOfModules) (i : M ≅ N) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj N ≅ 𝟙_ (PresheafModulesOn X V))
    (s : M.obj (op ⊤)) :
    trivValue N V e (i.hom.app (op ⊤) s)
      = trivValue M V ((restrictPresheafFunctor X V).mapIso i ≪≫ e) s := by
  simp only [trivValue, secOn_map M N i.hom V s]
  rfl

/-- ★★★★★**`X_s` は同型で移る**。

★機構は「自明化の集合が `e ↦ i|_V ≪≫ e` で 1 対 1 に対応する」ことだけである。 -/
theorem nonVanishing_of_iso (M N : X.PresheafOfModules) (i : M ≅ N) (s : M.obj (op ⊤)) :
    nonVanishing N (i.hom.app (op ⊤) s) = nonVanishing M s := by
  apply le_antisymm
  · intro x hx
    obtain ⟨V, e, hxe⟩ := (mem_nonVanishing_iff N _ x).1 hx
    rw [trivValue_of_iso M N i V e s] at hxe
    exact basicOpen_trivValue_le M V _ s hxe
  · intro x hx
    obtain ⟨V, e', hxe⟩ := (mem_nonVanishing_iff M s x).1 hx
    refine basicOpen_trivValue_le N V (((restrictPresheafFunctor X V).mapIso i).symm ≪≫ e') _ ?_
    rw [trivValue_of_iso M N i V _ s, basicOpen_trivValue_congr M V e' _ s]
    exact hxe

/-- ★★★★★★★★**空虚封じ** —— アフィンスキームの上では構造層 `𝒪_X` は ample である。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`n = 1`、`s ≝ λ⁻¹(1)`（`λ_ : 𝟙 ⊗ 𝟙 ≅ 𝟙` の逆で `1` を移したもの）と取ると
`X_s = X.basicOpen 1 = ⊤` で、`⊤` はアフィンである。
★★これで `IsAmple` が**満たされうる述語**であることが型で言える。 -/
theorem isAmple_unit [IsAffine X] : IsAmple (𝟙_ X.PresheafOfModules) := by
  refine ⟨inferInstance, fun x => ⟨1, one_pos, ?_⟩⟩
  refine ⟨(λ_ (𝟙_ X.PresheafOfModules)).inv.app (op ⊤) (1 : Γ(X, (⊤ : X.Opens))), ?_, ?_⟩
  · have h := nonVanishing_of_iso (𝟙_ X.PresheafOfModules)
      (𝟙_ X.PresheafOfModules ⊗ 𝟙_ X.PresheafOfModules)
      (λ_ (𝟙_ X.PresheafOfModules)).symm (1 : Γ(X, (⊤ : X.Opens)))
    show x ∈ nonVanishing (𝟙_ X.PresheafOfModules ⊗ 𝟙_ X.PresheafOfModules) _
    rw [show ((λ_ (𝟙_ X.PresheafOfModules)).symm).hom = (λ_ (𝟙_ X.PresheafOfModules)).inv from rfl]
      at h
    rw [h, nonVanishing_unit, Scheme.basicOpen_one]
    trivial
  · have h := nonVanishing_of_iso (𝟙_ X.PresheafOfModules)
      (𝟙_ X.PresheafOfModules ⊗ 𝟙_ X.PresheafOfModules)
      (λ_ (𝟙_ X.PresheafOfModules)).symm (1 : Γ(X, (⊤ : X.Opens)))
    show IsAffineOpen (nonVanishing (𝟙_ X.PresheafOfModules ⊗ 𝟙_ X.PresheafOfModules) _)
    rw [show ((λ_ (𝟙_ X.PresheafOfModules)).symm).hom = (λ_ (𝟙_ X.PresheafOfModules)).inv from rfl]
      at h
    rw [h, nonVanishing_unit, Scheme.basicOpen_one]
    exact isAffineOpen_top X

/-! ### ★出典の紐付け(`.src`)

★★`Proposition 1.4, (iv)` の証明が使う語彙である。
`ample ⟹ 射影埋め込み`（Serre の定理）は含まない。 -/

def nonVanishing.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(語彙——切断の非消失軌跡 X_s)",
    sectionId := "genell-prop-1-4" }

def basicOpen_trivValue_congr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(X_s は同じ V の上の自明化の取り方に依らない)",
    sectionId := "genell-prop-1-4" }

def presheafTensorPow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(語彙——前層加群のテンソル冪 M^{⊗n})",
    sectionId := "genell-prop-1-4" }

def IsAmple.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(語彙——ample の定義。Serre の定理は含まない)",
    sectionId := "genell-prop-1-4" }

def nonVanishing_of_iso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(X_s は同型で移ること)",
    sectionId := "genell-prop-1-4" }

def isAmple_unit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(空虚封じ——アフィンなら構造層は ample)",
    sectionId := "genell-prop-1-4" }

def sectionRatio.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(切断の比 s/t——同次座標の入口)",
    sectionId := "genell-prop-1-4" }

def sectionRatio_congr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(切断の比は自明化の取り方に依らないこと)",
    sectionId := "genell-prop-1-4" }

def nonVanishing.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "IsLocallyTrivial(可逆層の強い局所自明性)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.IsLocallyTrivial") 6,
    .citation "[ABC3]" "restrictPresheafFunctor(前層加群の開集合への制限)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.restrictPresheafFunctor") 6,
    .citation "[mathlib]" "AlgebraicGeometry.Scheme.basicOpen"
      (.inMathlib "AlgebraicGeometry.Scheme.basicOpen") 6,
    .otherPaper "[Stacks]" "28.27.1(ample の定義、tag 01PS)" 6,
    .implicitStep
      ("★mathlib に加群層の**茎**は 1 件も無い(2026-08-27 実測)ので、" ++
       "教科書の定義(茎で s_x ∉ 𝔪_x M_x)は書けない。" ++
       "★★代わりに局所自明化ごとに basicOpen を取って合併した" ++
       "——開であることが自動になり、自明化を選ばないので選択公理も要らない") 6,
    .implicitStep
      ("★★★異なる V, W にまたがる独立性" ++
       "(X_s ⊓ V = X.basicOpen (trivValue M V e s))は 2026-08-27 に閉じた" ++
       "——nonVanishing_inf。機構は trivValue_restrict(e.hom の自然性)だけである。" ++
       "★健全性も取れている——構造層なら X_s = Scheme.basicOpen s(nonVanishing_unit)。" ++
       "★★残るのは very ample(段 C2 の後)と Serre の定理(段 E)である") 6 ]

def nonVanishing_inf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(X_s は自明化の上でちょうど basicOpen であること)",
    sectionId := "genell-prop-1-4" }

def nonVanishing_unit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(健全性——構造層なら X_s = Scheme.basicOpen s)",
    sectionId := "genell-prop-1-4" }

def trivValue_restrict.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(自明化の制限と trivValue が可換であること)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
