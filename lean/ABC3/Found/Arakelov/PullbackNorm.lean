/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.PullbackMetric
import ABC3.Found.Arakelov.PullbackGen
import ABC3.Found.Arakelov.AMetricNorm
import ABC3.Meta.Claim

/-!
# 引き戻した切断のノルム —— `|f^*s|(q) = |s|(q ∘ f)`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

## ★★★★★★★★これは何か

`§9-742` で計量の引き戻し `AMetricPullback` が入り、`§9-741` で切断の引き戻し
`pullSec` が入った。本ファイルはその**両者を繋ぐ**:

    `|pullSec f s|_{f^*L̄}(q) = |s|_{L̄}(q ≫ f)`

★★これが**底変換のアルキメデス側**の核である——
`deg_K(L|_{Spec 𝓞_K}) = deg_F(L)` の証明で、
`Σ_τ log |pullSec s|_τ` を `Σ_σ log |s|_σ` に直すのに使う。

## ★★★★★機構（4 段、すべて在庫）

    pullback_h_pullTrivOfBase : 引き戻した計量の基準ノルムは元の基準ノルム（`§9-742`）
    pullSec_restrict          : 切断の引き戻しは制限と可換（`§9-741`）
    trivEquiv_pullSec         : 自明化での値は `f.c` で写る（`§9-741`）
    evalOn_pullback           : 関数の値は `q ≫ f` での値（在庫）

★`AMetric.norm_eq`（チャート独立性）で「`q ≫ f` の周りの `L` のチャート `W`」を選び、
その引き戻し `f⁻¹W` で測る、というのが全体の骨である。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
open ABC3.Found.GenEll

/-- ★★★★★★★★★★**引き戻した切断のノルムは、引き戻した点での元のノルムである**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

    `|pullSec f s|_{f^*L̄}(q) = |s|_{L̄}(q ≫ f)`

★★これが**底変換のアルキメデス側**の核である。 -/
theorem norm_pullSec {X Y : Scheme.{0}} (f : X ⟶ Y) (L : AMetric Y)
    (s : (L.sheaf.obj (op ⊤) : Type)) (q : Spec (CommRingCat.of ℂ) ⟶ X) :
    (AMetricPullback f L).norm (pullSec f L.sheaf ⊤ s) q = L.norm s (q ≫ f) := by
  obtain ⟨c⟩ := nonempty_normChart L.triv (q ≫ f)
  have hq : q ⁻¹ᵁ ((Opens.map f.base).obj c.V) = ⊤ := c.hp
  rw [AMetric.norm_eq (AMetricPullback f L) _ q ((Opens.map f.base).obj c.V)
       (pullTrivOfBase f L.sheaf c.V c.e le_rfl) hq,
      AMetric.norm_eq L s (q ≫ f) c.V c.e c.hp]
  show trivSecNorm ((pullbackPre f).obj L.sheaf) ((Opens.map f.base).obj c.V)
        (pullTrivOfBase f L.sheaf c.V c.e le_rfl) q hq (pullSec f L.sheaf ⊤ s)
      * (LocalMetric.pullback f L.triv L.metric).h ((Opens.map f.base).obj c.V)
          (pullTrivOfBase f L.sheaf c.V c.e le_rfl) q
    = trivSecNorm L.sheaf c.V c.e (q ≫ f) c.hp s * L.metric.h c.V c.e (q ≫ f)
  rw [pullback_h_pullTrivOfBase f L.triv L.metric c.V c.e le_rfl q hq]
  congr 1
  have hB : secOn ((pullbackPre f).obj L.sheaf) ((Opens.map f.base).obj c.V)
      (pullSec f L.sheaf ⊤ s) = pullSec f L.sheaf c.V (secOn L.sheaf c.V s) :=
    (pullSec_restrict f L.sheaf (le_top : c.V ≤ ⊤) s
      (le_top : (Opens.map f.base).obj c.V ≤ (Opens.map f.base).obj (⊤ : Y.Opens))).symm
  show ‖evalOn q _ hq (trivValue ((pullbackPre f).obj L.sheaf) _
      (pullTrivOfBase f L.sheaf c.V c.e le_rfl) (pullSec f L.sheaf ⊤ s))‖
    = ‖evalOn (q ≫ f) c.V c.hp (trivValue L.sheaf c.V c.e s)‖
  rw [trivValue_eq_trivEquiv, hB, trivEquiv_pullSec f L.sheaf c.V c.e (secOn L.sheaf c.V s)]
  congr 1

/-! ### ★出典の紐付け(`.src`) -/

def norm_pullSec.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(引き戻した切断のノルムは引き戻した点での元のノルム)",
    sectionId := "genell-def-1-1-ii" }

def norm_pullSec.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "pullback_h_pullTrivOfBase(引き戻した計量の基準ノルム、§9-742)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.pullback_h_pullTrivOfBase") 3,
    .citation "[ABC3]" "pullSec_restrict / trivEquiv_pullSec(切断の引き戻し、§9-741)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.trivEquiv_pullSec") 3,
    .citation "[ABC3]" "evalOn_pullback(関数の値は q ≫ f での値、在庫)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.evalOn_pullback") 3,
    .implicitStep
      ("★★これは**底変換のアルキメデス側**の核である。" ++
       "残るのは (a) 埋め込みの数え上げ(σ : F ↪ ℂ の延長はちょうど [K:F] 個)と、" ++
       "(b) **有限側**——Γ(f^*L) ≅ Γ(L) ⊗_{𝓞_F} 𝓞_K。" ++
       "★★★(b) は mathlib の PresheafOfModules.pullback が随伴の左随伴として" ++
       "**抽象的にしか定義されていない**(各開集合ごとの式が無い)ため、" ++
       "現状の設計では直接には書けない。2026-08-28 実測") 4 ]

end ABC3.Found.Arakelov
