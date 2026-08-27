/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.AmpleDef
import ABC3.Meta.Claim

/-!
# 自明化の値は**テンソル積で掛け算になる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★台帳 `arakelov-coherent-base-metric` の**鍵の恒等式**

2026-08-28 に `Definition 1.1` の項目全体の `.src` を下げた。理由は

    `TorsorMetric.base` は対象ごとの `Classical.choice` なので
    `base_{[L·M]} ≠ base_{[L]} ⊗ base_{[M]}` となり、
    群法則が**計量のテンソル積を表さない**

であった。★塡がりの本体は「**計量が掛け算になる**」を言えるようにすることである。

## ★★★★★★本ファイルが与えるもの

    `trivValue (A ⊗ B) V (tensorTriv eA eB) (s ⊗ₜ t)`
      `= trivValue A V eA s * trivValue B V eB t`

★★★**これは `rfl` である**（2026-08-28 実測）。機構は 3 つとも mathlib で `rfl`:

| 段 | 根拠 |
|---|---|
| 制限のモノイダル構造 `μIso` | ★`Iso.refl`（`PushforwardZeroMonoidal.lean`） |
| `(eA ⊗ᵢ eB)` の切断での値 | ★★`tensorHom_app`（`rfl`） |
| `λ_ (𝟙_)` の切断での値 | ★★★`leftUnitor_hom_app`（`rfl`）——`r ⊗ₜ m ↦ r • m` |

★`set_option backward.isDefEq.respectTransparency false` が要る
（`tools\lean-idioms.md`——mathlib 自身が同じ設定を使っている）。

## ★★★★★なぜこれが鍵か

局所自明な前層加群の計量は、**自明化ごとの基準ノルム** `h_{V,e} ≝ ‖e⁻¹(1)‖` で表せる。
そこでは切断のノルムが

    `‖s‖(p) = |trivValue F V e s (p)| · h_{V,e}(p)`

と書け、★★テンソル積の計量は `h^{A⊗B}_{V, eA⊗eB} = h^A_{V,eA} · h^B_{V,eB}` で入る。
★★★本ファイルの恒等式が「そのとき切断のノルムも掛け算になる」を保証する
——それが `Classical.choice` の基準計量では**落ちていた**加法性である。

## ★残っている段（明示）

★★**ファイバーの段の橋**——`arcFiber p (A ⊗ B) ≅ arcFiber p A ⊗_ℂ arcFiber p B`——は
**まだ無い**。mathlib の `SheafOfModules` にモノイダル構造が無く（2026-08-28 実測）、
引き戻しのモノイダル性も無いからである。
★本ファイルは**前層の段**で恒等式を取り、その橋を迂回する道を開く。

## ★★ample の側への副産物

    `X.basicOpen (trivValue (A ⊗ B) V (tensorTriv eA eB) (s ⊗ₜ t))`
      `= X.basicOpen (trivValue A V eA s) ⊓ X.basicOpen (trivValue B V eB t)`

すなわち **`X_{s⊗t} = X_s ∩ X_t`**。★これは very ample の議論
（`n·D` の切断を `s^{⊗n}` で作る段）で直に要る。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite
open ABC3.Found.GenEll

variable {X : Scheme.{u}}

/-! ## ★★★テンソル積の自明化 -/

/-- ★★**テンソル積の自明化**——`IsLocallyTrivial.tensor` が使っているものと同じ。

★`(A ⊗ B)|_V ≅ A|_V ⊗ B|_V ≅ 𝟙_ ⊗ 𝟙_ ≅ 𝟙_`。 -/
noncomputable def tensorTriv {A B : X.PresheafOfModules} {V : X.Opens}
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V)) :
    (restrictPresheafFunctor X V).obj (A ⊗ B) ≅ 𝟙_ (PresheafModulesOn X V) :=
  (restrictPresheafTensor A B).symm ≪≫ (eA ⊗ᵢ eB) ≪≫ λ_ (𝟙_ (PresheafModulesOn X V))

/-- ★**制限は `⊗ₜ` を `⊗ₜ` へ送る**。

★★前層のテンソル積は切断ごとに `Γ(A,W) ⊗_{Γ(X,W)} Γ(B,W)` **そのもの**なので
（mathlib の `tensorObj` の `obj X := M₁.obj X ⊗ M₂.obj X`）、これは `rfl` である。 -/
theorem secOn_tmul {A B : X.PresheafOfModules} (V : X.Opens)
    (s : A.obj (op ⊤)) (t : B.obj (op ⊤)) :
    secOn (A ⊗ B) V (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type u)] t)
      = secOn A V s ⊗ₜ[(Γ(X, V) : Type u)] secOn B V t := rfl

/-! ## ★★★★★★★★★鍵の恒等式 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★**自明化の値はテンソル積で掛け算になる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★**`rfl` である**（2026-08-28 実測）——3 つの機構がすべて mathlib で `rfl` だから:
制限の `μIso` は `Iso.refl`、`tensorHom_app` は `rfl`、
`leftUnitor_hom_app` は `r ⊗ₜ m ↦ r • m` で `rfl`。

★★★これが台帳 `arakelov-coherent-base-metric` の**鍵**である
——「計量が掛け算になる」を前層の段で言えるようにする。 -/
theorem trivValue_tensor {A B : X.PresheafOfModules} {V : X.Opens}
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V))
    (s : A.obj (op ⊤)) (t : B.obj (op ⊤)) :
    trivValue (A ⊗ B) V (tensorTriv eA eB) (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type u)] t)
      = trivValue A V eA s * trivValue B V eB t := rfl

/-! ## ★★ample の側への副産物 -/

/-- ★★★★**`X_{s⊗t} = X_s ∩ X_t`**（自明化を固定した形）。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★very ample の議論（`n·D` の切断を `s^{⊗n}` で作る段）で直に要る。 -/
theorem basicOpen_trivValue_tensor {A B : X.PresheafOfModules} {V : X.Opens}
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V))
    (s : A.obj (op ⊤)) (t : B.obj (op ⊤)) :
    X.basicOpen (trivValue (A ⊗ B) V (tensorTriv eA eB)
        (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type u)] t))
      = X.basicOpen (trivValue A V eA s) ⊓ X.basicOpen (trivValue B V eB t) := by
  rw [trivValue_tensor, Scheme.basicOpen_mul]

/-- ★★★**非消失軌跡の側**——各自明化の対が与える交わりは `X_{s⊗t}` に含まれる。 -/
theorem le_nonVanishing_tensor {A B : X.PresheafOfModules} {V : X.Opens}
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V))
    (s : A.obj (op ⊤)) (t : B.obj (op ⊤)) :
    X.basicOpen (trivValue A V eA s) ⊓ X.basicOpen (trivValue B V eB t)
      ≤ nonVanishing (A ⊗ B) (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type u)] t) := by
  rw [← basicOpen_trivValue_tensor eA eB s t]
  exact basicOpen_trivValue_le _ V (tensorTriv eA eB) _

/-! ### ★出典の紐付け(`.src`) -/

def tensorTriv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——テンソル積の自明化)",
    sectionId := "genell-def-1-1-i" }

def trivValue_tensor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——自明化の値がテンソル積で掛け算になること)",
    sectionId := "genell-def-1-1-i" }

def basicOpen_trivValue_tensor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(X_{s⊗t} = X_s ∩ X_t)",
    sectionId := "genell-prop-1-4" }

def trivValue_tensor.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "PresheafOfModules.Monoidal.tensorHom_app / leftUnitor_hom_app(ともに rfl)"
      (.inMathlib "PresheafOfModules.Monoidal.tensorHom_app") 3,
    .citation "[mathlib]" "pushforward₀OfCommRingCat の Monoidal(μIso = Iso.refl)"
      (.inMathlib "PresheafOfModules.pushforward₀OfCommRingCat") 3,
    .citation "[ABC3]" "trivValue(自明化のもとで切断が定める関数)"
      (.inProject "ABC3" "ABC3.Found.GenEll.trivValue") 3,
    .implicitStep
      ("★★★★★これは台帳 arakelov-coherent-base-metric の**鍵の恒等式**である。" ++
       "2026-08-28 に Definition 1.1 の項目全体の .src を下げた理由は " ++
       "TorsorMetric.base が Classical.choice で base_{[L·M]} ≠ base_{[L]} ⊗ base_{[M]} " ++
       "となることであり、塡がりの本体は「計量が掛け算になる」を言えるようにすることであった") 3,
    .implicitStep
      ("★★残っているのは**ファイバーの段の橋** " ++
       "arcFiber p (A ⊗ B) ≅ arcFiber p A ⊗_ℂ arcFiber p B である。" ++
       "mathlib の SheafOfModules にモノイダル構造が無く(2026-08-28 実測)、" ++
       "引き戻しのモノイダル性も無い。★本ファイルは前層の段で恒等式を取り、その橋を迂回する道を開く") 3 ]

end ABC3.Found.Arakelov
