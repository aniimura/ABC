/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.AMetricHom
import ABC3.Found.Arakelov.PicLTPull
import ABC3.Meta.Claim

/-!
# 引き戻した層の**自明化**を名前にする（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

## ★★★★★★★`Definition 1.1, (ii)` —— 台の側は**すべて在庫**である

2026-08-28 に台帳 `arakelov-pullback-monoidal` を立てたのは**在庫確認の怠り**であった。
★以下はすべて `Found/Arakelov/` に `sorry` 無しで在る:

| 在庫 | 内容 |
|---|---|
| `pullbackPre f` | スキームの射に沿った**前層**加群の引き戻し |
| `pullbackPreOplax` | ★`(pullbackPre f).OplaxMonoidal` |
| `pullbackDelta` / `pullbackEta` | `f^*(P⊗Q) ⟶ f^*P ⊗ f^*Q` と `f^*𝒪 ⟶ 𝒪` |
| `isLocallyTrivial_pullbackPre` | ★★局所自明性は引き戻しで保たれる |
| `bcIso` | ★★★Beck–Chevalley `(f|)^*(L|_W) ≅ (f^*L)|_{f⁻¹W}` |
| `pullbackOnUnitIso` | `(f|)^* 𝟙_ ≅ 𝟙_` |
| `pullbackUnitIso` / `pullbackSheafTensorIso'` | 層の段の単位・テンソル |

## ★★★★★★本ファイルがすること

`isLocallyTrivial_pullbackPre` の証明の中に**自明化の輸送が calc として書かれている**が、
名前が付いていないので**外から使えない**。★本ファイルはそれを `pullTrivOfBase` として
取り出す——計量を運ぶ段（次のブロック）はこれを名前で呼ぶ。

## ★★★★★★★★★次のブロックの設計（明示）

`AMetric` の引き戻しは `LocalMetric.tensor`（`TensorMetric.lean`）と**同じ型**である:

1. ★**チャート**: `V : X.Opens`・`p` に対し、`L|_W ≅ 𝟙_` となる `W : Y.Opens` と
   `V' ≤ V ⊓ f⁻¹W` で `p ⁻¹ᵁ V' = ⊤` なるものを取る
   （存在は `L.triv` ＋ `transportSieve`）。
2. ★★**候補値**:

       `‖evalOn p V' (transUnit (f^*L) V' (pullTrivOfBase …) (e|_{V'}))‖⁻¹`
         `· L.metric.h W eW (p ≫ f)`

3. ★★★**チャート独立性**の核は
   「`transUnit` が輸送 `pullTrivOfBase` と可換」——すなわち
   `transUnit (f^*L) V' (pullTrivOfBase eW₁) (pullTrivOfBase eW₂)`
   が `transUnit L W eW₁ eW₂` の**引き戻し**であること。
   ★これは `LocalMetric.tensor` の `transUnit_tensorTriv` に当たる補題である。
4. ★★★★あとは `chartH_triv_indep` / `chartH_shrink` / `chartH_indep` と
   同じ 3 段を書けばよい。

★★★★★見積もり **2〜4 ブロック**（`LocalMetric.tensor` と同規模）。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory TopologicalSpace MonoidalCategory Opposite

/-- ★★★★★★**引き戻した層の自明化**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

`L|_W ≅ 𝟙_`（`Y` の上）から `(f^*L)|_{V'} ≅ 𝟙_`（`X` の上、`V' ≤ f⁻¹W`）を作る。

★機構は `isLocallyTrivial_pullbackPre` の証明にある 4 段そのもの:
Beck–Chevalley（`bcIso`）→ `eW` の輸送 → `(f|)^*𝟙_ ≅ 𝟙_`（`pullbackOnUnitIso`）→
制限の推移律。★★本ファイルはそれに**名前を付ける**——
計量を運ぶ段が外から呼べるようにするためである。 -/
noncomputable def pullTrivOfBase {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (W : Y.Opens) (eW : (restrictPresheafFunctor Y W).obj L ≅ 𝟙_ (PresheafModulesOn Y W))
    {V' : X.Opens} (hV'W : V' ≤ (Opens.map f.base).obj W) :
    (restrictPresheafFunctor X V').obj ((pullbackPre f).obj L)
      ≅ 𝟙_ (PresheafModulesOn X V') :=
  (Iso.refl _ : (restrictPresheafFunctor X V').obj ((pullbackPre f).obj L)
      ≅ (restrictOnFunctor hV'W).obj
          ((restrictPresheafFunctor X ((Opens.map f.base).obj W)).obj ((pullbackPre f).obj L)))
    ≪≫ (restrictOnFunctor hV'W).mapIso (bcIso f W L).symm
    ≪≫ (restrictOnFunctor hV'W).mapIso ((pullbackPreOn f W).mapIso eW)
    ≪≫ (restrictOnFunctor hV'W).mapIso (pullbackOnUnitIso f W)
    ≪≫ (Iso.refl _ : (restrictOnFunctor hV'W).obj
          (𝟙_ (PresheafModulesOn X ((Opens.map f.base).obj W))) ≅ 𝟙_ (PresheafModulesOn X V'))

/-! ### ★出典の紐付け(`.src`) -/

def pullTrivOfBase.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(引き戻した層の自明化——輸送を名前にしたもの)",
    sectionId := "genell-def-1-1-ii" }

def pullTrivOfBase.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "bcIso(Beck–Chevalley)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.bcIso") 3,
    .citation "[ABC3]" "pullbackOnUnitIso((f|)^* 𝟙_ ≅ 𝟙_)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.pullbackOnUnitIso") 3,
    .citation "[ABC3]" "isLocallyTrivial_pullbackPre(局所自明性は引き戻しで保たれる)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.isLocallyTrivial_pullbackPre") 3,
    .implicitStep
      ("★★★★★**訂正(2026-08-28)**: 台帳 arakelov-pullback-monoidal を立てたのは" ++
       "在庫確認の怠りであった。mathlib に Monoidal インスタンスが無いのは事実だが、" ++
       "本プロジェクトは pullbackPreOplax を既に建てており、" ++
       "自明化の輸送も isLocallyTrivial_pullbackPre の証明の中にある") 3,
    .implicitStep
      ("★★次のブロックの設計: AMetric の引き戻しは LocalMetric.tensor と同じ型である。" ++
       "(1) チャート(L|_W ≅ 𝟙_ となる W と V' ≤ V ⊓ f⁻¹W)、" ++
       "(2) 候補値 ‖evalOn(transUnit (f^*L) V' (pullTrivOfBase …) (e|_{V'}))‖⁻¹ · h W eW (p ≫ f)、" ++
       "(3) チャート独立性の核は「transUnit が輸送と可換」" ++
       "(LocalMetric.tensor の transUnit_tensorTriv に当たる補題)。" ++
       "★見積もり 2〜4 ブロック") 3 ]

end ABC3.Found.Arakelov
