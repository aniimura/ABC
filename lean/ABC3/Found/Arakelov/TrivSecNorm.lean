/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.TrivTensor
import ABC3.Found.Arakelov.ArcSemilinear
import ABC3.Meta.Claim

/-!
# 切断のノルムは**テンソル積で掛け算になる**——ファイバーを経由しない道（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★台帳 `arakelov-coherent-base-metric` の段 2

段 1（`TrivTensor.lean`）で `trivValue` がテンソル積で掛け算になることを取った。
★本ファイルはそれを**切断のノルム**に載せる。

## ★★★★★★★★ファイバーの橋を**迂回する**

`arcFiber p (A ⊗ B) ≅ arcFiber p A ⊗_ℂ arcFiber p B` は mathlib に無く
（`SheafOfModules` にモノイダル構造が無い、2026-08-28 実測）、作るのは PR 規模である。

★★**しかし切断のノルムはファイバーを経由しなくても書ける**:

    `‖s‖_{V,e}(p) ≝ ‖evalOn p V hp (trivValue F V e s)‖`

これは「自明化 `e` の生成切断のノルムを 1 に正規化したときの `s` のノルム」であり、
★★★在庫の `genNorm_arcEvalOnTop` が

    `genNorm F hF V g p h (arcEvalOnTop F p V h (c • g)) = ‖evalOn p V h c‖`

としてすでに保証している——**ファイバーで測っても同じ値になる**。

## ★★★★★★★★★そして掛け算になる

    `‖s ⊗ t‖_{V, eA⊗eB} = ‖s‖_{V,eA} · ‖t‖_{V,eB}`

★機構は 2 つだけ: `trivValue_tensor`（段 1）と **`evalOn` が環準同型であること**。

★★★これが「`Classical.choice` の基準計量では落ちていた加法性」の**本体**である。
一般の計量は基準ノルム `h_{V,e} > 0` を掛けた

    `‖s‖(p) = ‖evalOn p V hp (trivValue F V e s)‖ · h_{V,e}(p)`

であり、テンソル積は `h_{A⊗B, eA⊗eB} ≝ h_{A,eA} · h_{B,eB}` で入る。

## ★残っている段（明示）

1. ★`h_{V,e}` を担ぐ構造体（自明化の取り替えで `|u|` で割られる整合条件つき）
2. ★★その整合条件のもとで `‖s‖(p)` が `(V, e)` に依らないこと
3. ★★★連続性と、`APic` の群法則への載せ替え
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite
open ABC3.Found.GenEll

variable {X : Scheme.{0}}

/-! ## ★★`evalOn` は環準同型である -/

/-- ★`evalOn` は積を保つ。 -/
theorem evalOn_mul (p : Spec (CommRingCat.of ℂ) ⟶ X) (V : X.Opens) (hp : p ⁻¹ᵁ V = ⊤)
    (c c' : ((X.presheaf.obj (op V)) : Type)) :
    evalOn p V hp (c * c') = evalOn p V hp c * evalOn p V hp c' := by
  show (Scheme.ΓSpecIso (CommRingCat.of ℂ)).hom.hom
      (((Spec (CommRingCat.of ℂ)).presheaf.map (homOfLE (le_of_eq hp.symm)).op).hom
        ((p.app V).hom (c * c'))) = _
  rw [map_mul, map_mul, map_mul]
  rfl

/-! ★`evalOn p V hp 1 = 1` は在庫（`ArcLocalMetric.lean` の `evalOn_one`）。 -/

/-! ## ★★★★★★自明化つきの切断ノルム -/

/-- ★★★★**自明化 `e` で正規化した切断のノルム**——ファイバーを経由しない。

★在庫の `genNorm_arcEvalOnTop` が「ファイバーで測っても同じ値」を保証している。 -/
noncomputable def trivSecNorm (F : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) (s : F.obj (op ⊤)) : ℝ :=
  ‖evalOn p V hp (trivValue F V e s)‖

theorem trivSecNorm_nonneg (F : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) (s : F.obj (op ⊤)) :
    0 ≤ trivSecNorm F V e p hp s :=
  norm_nonneg _

/-- ★自明化を取り替えると単元の絶対値が掛かる。 -/
theorem trivSecNorm_congr (F : X.PresheafOfModules) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤) :
    ∃ u : ((X.presheaf.obj (op V)) : Type), IsUnit u ∧
      ∀ s : F.obj (op ⊤),
        trivSecNorm F V e' p hp s = trivSecNorm F V e p hp s * ‖evalOn p V hp u‖ := by
  obtain ⟨u, hu, hus⟩ := trivValue_congr' F V e e'
  refine ⟨u, hu, fun s => ?_⟩
  show ‖evalOn p V hp (trivValue F V e' s)‖ = _
  rw [hus s, evalOn_mul, norm_mul]
  rfl

/-! ## ★★★★★★★★★掛け算になること -/

/-- ★★★★★★★★★**切断のノルムはテンソル積で掛け算になる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★機構は 2 つだけ——`trivValue_tensor`（段 1）と `evalOn` が環準同型であること。

★★★これが台帳 `arakelov-coherent-base-metric` の本体である
——`Classical.choice` の基準計量では落ちていた加法性が、
自明化の側で**恒等式として**出る。 -/
theorem trivSecNorm_tensor {A B : X.PresheafOfModules} (V : X.Opens)
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤)
    (s : A.obj (op ⊤)) (t : B.obj (op ⊤)) :
    trivSecNorm (A ⊗ B) V (tensorTriv eA eB) p hp
        (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t)
      = trivSecNorm A V eA p hp s * trivSecNorm B V eB p hp t := by
  show ‖evalOn p V hp (trivValue (A ⊗ B) V (tensorTriv eA eB) _)‖ = _
  rw [trivValue_tensor, evalOn_mul, norm_mul]
  rfl

/-- ★★**ノルムの対数は足し算になる**——Green 関数の側の形。

★どちらの切断も `p` で消えていないときである（`log 0` を避ける）。 -/
theorem log_trivSecNorm_tensor {A B : X.PresheafOfModules} (V : X.Opens)
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V))
    (eB : (restrictPresheafFunctor X V).obj B ≅ 𝟙_ (PresheafModulesOn X V))
    (p : Spec (CommRingCat.of ℂ) ⟶ X) (hp : p ⁻¹ᵁ V = ⊤)
    (s : A.obj (op ⊤)) (t : B.obj (op ⊤))
    (hs : trivSecNorm A V eA p hp s ≠ 0) (ht : trivSecNorm B V eB p hp t ≠ 0) :
    Real.log (trivSecNorm (A ⊗ B) V (tensorTriv eA eB) p hp
        (s ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] t))
      = Real.log (trivSecNorm A V eA p hp s) + Real.log (trivSecNorm B V eB p hp t) := by
  rw [trivSecNorm_tensor, Real.log_mul hs ht]

/-! ### ★出典の紐付け(`.src`) -/

def evalOn_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——evalOn は環準同型)",
    sectionId := "genell-def-1-1-i" }

def trivSecNorm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(自明化で正規化した切断のノルム。ファイバーを経由しない)",
    sectionId := "genell-def-1-1-i" }

def trivSecNorm_tensor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(切断のノルムがテンソル積で掛け算になること)",
    sectionId := "genell-def-1-1-i" }

def trivSecNorm_tensor.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "trivValue_tensor(自明化の値がテンソル積で掛け算になる)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.trivValue_tensor") 3,
    .citation "[ABC3]" "genNorm_arcEvalOnTop(正規化ノルムは正則関数の絶対値になる)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.genNorm_arcEvalOnTop") 3,
    .implicitStep
      ("★★★★★ファイバーの橋 arcFiber p (A ⊗ B) ≅ arcFiber p A ⊗_ℂ arcFiber p B は " ++
       "mathlib に無く(SheafOfModules にモノイダル構造が無い、2026-08-28 実測)作るのは PR 規模である。" ++
       "★本ファイルは切断のノルムを自明化の側で書くことで**その橋を迂回した**。" ++
       "在庫の genNorm_arcEvalOnTop がファイバーで測っても同じ値になることを保証している") 3,
    .implicitStep
      ("★★残る段: (1) 基準ノルム h_{V,e} を担ぐ構造体(自明化の取り替えで |u| で割られる整合条件つき)、" ++
       "(2) その整合条件のもとで ‖s‖(p) が (V,e) に依らないこと、" ++
       "(3) 連続性と APic の群法則への載せ替え") 3 ]

end ABC3.Found.Arakelov
