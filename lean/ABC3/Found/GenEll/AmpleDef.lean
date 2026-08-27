/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.PicLocalTrivial
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

## ★残っている段（明示）

★★**異なる `V, W` にまたがる独立性**——`X_s ⊓ V = X.basicOpen (trivValue M V e s)`——
はまだ無い。要るのは「自明化の制限と `trivValue` が可換」であること
（`Found/Arakelov/PicBasicTrivial.lean` の `trivialOfLe` の自然性）である。
★★★`ample` の基本性質（`X_s` が `L^{⊗n}` の取り方に依らない等）はその先である。
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
      ("★★★残っているのは**異なる V, W にまたがる独立性**" ++
       "(X_s ⊓ V = X.basicOpen (trivValue M V e s))である。" ++
       "要るのは自明化の制限(trivialOfLe)と trivValue が可換であること") 6 ]

end ABC3.Found.GenEll
