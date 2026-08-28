/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.SpanPullSec
import ABC3.Found.Arakelov.TensorSurj
import ABC3.Found.Arakelov.InvertibleIndex
import ABC3.Found.Arakelov.PullbackTensor
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★`Γ(f^*L)` は**引き戻した大域切断で生成される**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

## ★★★★★★★★★★これは何か——台帳 `presheaf-pullback-global-sections` を塞ぐ

    `Γ(f^*L) = Γ(X,⊤)-span { pullSec f L ⊤ s : s ∈ Γ(L) }`   （`L` は算術直線束、`f` は任意）

★mathlib の `PresheafOfModules.pullback` は**随伴の左随伴としてしか定義されていない**
（各開集合ごとの式が無い）ので、`Γ(f^*L)` は直接には計算できない。
★★本ファイルは**テンソルで回して**その計算を回避する。

## ★★★★★機構（3 段、材料はすべて在庫）

    span_pullSecT_unit  : 単位束では成り立つ（`§9-789`）
    span_pullSecT_congr : 同型で移せる（`§9-789`）
    ★ここから `A ⊗ A⁻¹ ≅ Ō` で `A ⊗ B` の場合に持ち上げ、
    delta_pullSec       : `δ(pullSec (a ⊗ b)) = pullSec a ⊗ pullSec b`（`§9-743`）
    isIso_pullbackDelta : `δ` は同型（在庫）
    surjective_of_map_surjective : `α ⊗ β` 全射 ⟹ `α` 全射（`§9-785`）

★★★筋: `span (range (pullSec_{A⊗B})) = ⊤` を `δ` で送ると
`span { pullSec a ⊗ pullSec b } = ⊤` になり、
これは `P ⊗ Q → Γ(f^*A) ⊗ Γ(f^*B)`（`P`, `Q` は span）が**全射**であることを言う。
★★★★`§9-785` でその因子を落とすと `P = ⊤`、すなわち求める生成性が出る。

## ★測定の記録

**`f` の全射性は要らない**——単位束の場合が `pullbackUnitPreIso`（任意の `f`）で片付き、
テンソルの議論も `f` に条件を課さないからである。
★最初の見立て（左 Kan 拡大の colimit が潰れるには `f` 全射が要る）は**不要な仮定**であった。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
open ABC3.Found.GenEll
open scoped TensorProduct

/-! ## ★切断の引き戻しの加法性 -/

theorem pullSecT_add' {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules)
    (s t : (L.obj (op ⊤) : Type)) :
    pullSecT f L (s + t) = pullSecT f L s + pullSecT f L t :=
  map_add (((PresheafOfModules.pullbackPushforwardAdjunction
    (pullbackPhi f)).unit.app L).app (op ⊤)).hom s t

theorem pullSecT_zero' {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Y.PresheafOfModules) :
    pullSecT f L 0 = 0 :=
  map_zero (((PresheafOfModules.pullbackPushforwardAdjunction
    (pullbackPhi f)).unit.app L).app (op ⊤)).hom

/-! ## ★★★★★★★★★★テンソルの因子を落とす -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★**`A ⊗ B` で生成性が成り立てば `A` でも成り立つ**（`B` の側が可逆なら）。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★機構: 在庫の `delta_pullSec`（`§9-743`）で
`δ(pullSec (a ⊗ b)) = pullSec a ⊗ pullSec b`。
`δ` は同型（`isIso_pullbackDelta`）なので `span` は `⊤` のまま移り、
`P ⊗ Q → Γ(f^*A) ⊗ Γ(f^*B)` が全射になる。
★★あとは `§9-785` の `surjective_of_map_surjective` で因子を落とす。 -/
theorem span_pullSecT_of_mul {X Y : Scheme.{0}} (f : X ⟶ Y) (A B : Y.PresheafOfModules)
    [Module.Invertible (Γ(X, (⊤ : X.Opens)) : Type)
      ((((pullbackPre f).obj B).obj (op (⊤ : X.Opens))) : Type)]
    (hAB : Submodule.span (Γ(X, (⊤ : X.Opens)) : Type)
      (Set.range (pullSecT f (A ⊗ B))) = ⊤) :
    Submodule.span (Γ(X, (⊤ : X.Opens)) : Type) (Set.range (pullSecT f A)) = ⊤ := by
  haveI := isIso_pullbackDelta f A B
  set δlin : ((((pullbackPre f).obj (A ⊗ B)).obj (op (⊤ : X.Opens))) : Type)
      →ₗ[(Γ(X, (⊤ : X.Opens)) : Type)]
      ((((pullbackPre f).obj A ⊗ (pullbackPre f).obj B).obj (op (⊤ : X.Opens))) : Type) :=
    (((pullbackDelta f A B).app (op (⊤ : X.Opens))).hom) with hδ
  have hδsurj : Function.Surjective δlin :=
    ((PresheafOfModules.evaluation X.ringCatSheaf.obj (op (⊤ : X.Opens))).mapIso
      (asIso (pullbackDelta f A B))).toLinearEquiv.surjective
  set P := Submodule.span (Γ(X, (⊤ : X.Opens)) : Type) (Set.range (pullSecT f A)) with hP
  set Q := Submodule.span (Γ(X, (⊤ : X.Opens)) : Type) (Set.range (pullSecT f B)) with hQ
  set m := TensorProduct.map (P.subtype) (Q.subtype) with hm
  have hmem : ∀ x : ((A.obj (op (⊤ : Y.Opens)) : Type)
      ⊗[(Γ(Y, (⊤ : Y.Opens)) : Type)] (B.obj (op (⊤ : Y.Opens)) : Type)),
      δlin (pullSecT f (A ⊗ B) x) ∈ LinearMap.range m := by
    intro x
    induction x using TensorProduct.induction_on with
    | zero =>
        have h0 : pullSecT f (A ⊗ B) (0 : ((A.obj (op (⊤ : Y.Opens)) : Type)
            ⊗[(Γ(Y, (⊤ : Y.Opens)) : Type)] (B.obj (op (⊤ : Y.Opens)) : Type))) = 0 :=
          pullSecT_zero' f (A ⊗ B)
        rw [h0, map_zero]
        exact Submodule.zero_mem _
    | tmul a b =>
        refine ⟨(⟨pullSecT f A a, Submodule.subset_span ⟨a, rfl⟩⟩ : P)
          ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)]
          (⟨pullSecT f B b, Submodule.subset_span ⟨b, rfl⟩⟩ : Q), ?_⟩
        show (pullSecT f A a) ⊗ₜ[(Γ(X, (⊤ : X.Opens)) : Type)] (pullSecT f B b) = _
        exact (delta_pullSec f A B ⊤ a b).symm
    | add u v hu hv =>
        have ha : pullSecT f (A ⊗ B) (u + v)
            = pullSecT f (A ⊗ B) u + pullSecT f (A ⊗ B) v := pullSecT_add' f (A ⊗ B) u v
        rw [ha, map_add]
        exact Submodule.add_mem _ hu hv
  have hspan2 : Submodule.span (Γ(X, (⊤ : X.Opens)) : Type)
      ((δlin : _ → _) '' (Set.range (pullSecT f (A ⊗ B)))) = ⊤ := by
    rw [← Submodule.map_span, hAB, Submodule.map_top, LinearMap.range_eq_top]
    exact hδsurj
  have hle : Submodule.span (Γ(X, (⊤ : X.Opens)) : Type)
      ((δlin : _ → _) '' (Set.range (pullSecT f (A ⊗ B)))) ≤ LinearMap.range m := by
    rw [Submodule.span_le]
    rintro z ⟨w, ⟨x, rfl⟩, rfl⟩
    exact hmem x
  have htop : LinearMap.range m = ⊤ := eq_top_iff.mpr (hspan2 ▸ hle)
  have hsurjm : Function.Surjective m := LinearMap.range_eq_top.mp htop
  have hsurjP : Function.Surjective (P.subtype) :=
    surjective_of_map_surjective (Γ(X, (⊤ : X.Opens)) : Type) P _ Q _ P.subtype Q.subtype hsurjm
  have hres := LinearMap.range_eq_top.mpr hsurjP
  rwa [Submodule.range_subtype] at hres

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★**`Γ(f^*L)` は引き戻した大域切断で生成される**。

原文 (GenEll p.3):
> (ii) Let φ : Y → X be a morphism of normal, Z-proper, Z-flat schemes. Then there is an evident notion of pull-back of arithmetic line bundles by φ.

★★★これが台帳 `presheaf-pullback-global-sections` の**閉じ方**である
——`Γ(f^*L)` を計算せずに、`A ⊗ A⁻¹ ≅ Ō` と単位束の場合で回した。 -/
theorem span_pullSecT_of_inv {X Y : Scheme.{0}} (f : X ⟶ Y) (L : AInv Y) :
    Submodule.span (Γ(X, (⊤ : X.Opens)) : Type)
      (Set.range (pullSecT f L.carrier.sheaf)) = ⊤ := by
  obtain ⟨φ, -⟩ := L.isInv
  haveI hinvB : Module.Invertible (Γ(X, (⊤ : X.Opens)) : Type)
      ((((pullbackPre f).obj L.inv.sheaf).obj (op (⊤ : X.Opens))) : Type) :=
    invertible_gammaPre (AInv.pullback f L.symm)
  refine span_pullSecT_of_mul f L.carrier.sheaf L.inv.sheaf ?_
  exact span_pullSecT_congr f φ.symm (span_pullSecT_unit f)

/-! ### ★出典の紐付け(`.src`) -/

def span_pullSecT_of_inv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (ii)(Γ(f^*L) は引き戻した大域切断で生成される)",
    sectionId := "genell-def-1-1-ii" }

def span_pullSecT_of_inv.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "delta_pullSec(δ は純テンソルを純テンソルに送る、§9-743)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.delta_pullSec") 3,
    .citation "[ABC3]" "isIso_pullbackDelta(δ は同型、在庫)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.isIso_pullbackDelta") 3,
    .citation "[ABC3]" "span_pullSecT_unit / span_pullSecT_congr(§9-789)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.span_pullSecT_unit") 3,
    .citation "[ABC3]" "surjective_of_map_surjective(§9-785)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.surjective_of_map_surjective") 4,
    .implicitStep
      ("★★★これで台帳 presheaf-pullback-global-sections が塞がった。" ++
       "mathlib の PresheafOfModules.pullback は随伴の左随伴としてしか定義されていないが、" ++
       "Γ(f^*L) を**計算せずに** A ⊗ A⁻¹ ≅ Ō と単位束の場合で回した") 4,
    .implicitStep
      ("★測定: f の全射性は要らない。単位束の場合が pullbackUnitPreIso(任意の f)で片付き、" ++
       "テンソルの議論も f に条件を課さないからである") 3 ]

end ABC3.Found.Arakelov
