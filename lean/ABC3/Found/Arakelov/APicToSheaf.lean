/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.SheafifyTriv
import ABC3.Found.Arakelov.PicEvalIso
import ABC3.Found.Arakelov.PicSheafGroup
import ABC3.Found.Arakelov.PicClassGroup
import ABC3.Found.Arakelov.PicAssoc
import ABC3.Found.Arakelov.PicGammaInv
import ABC3.Meta.Claim

/-!
# **`APic`（前層・計量つき）から `Pic`（層）へ**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★★★台帳 `arakelov-degF-finite-places` の**段 A3**

無条件の `deg_F`（[Szp] を仮定しない形）に要るのは

    `ADiv(F)/APrc(F) ≅ APic(Spec 𝓞_F)`

だけになった（`§9-746`・`§9-747`）。その**有限素点側**は

    `APicM (Spec 𝓞_F) → PicSheaf (Spec 𝓞_F) ≃* Cl(F)`

で読める。★★本ファイルはその**第 1 の矢印**を作る。

| 段 | 内容 | 状態 |
|---|---|---|
| A1 | 層化は局所自明性を保つ | ✅ 在庫 |
| A2 | 自明化・遷移単元の輸送 | ✅ `SheafifyTriv.lean`（§9-748） |
| A3 | **`APicM X →* PicSheaf X`** | ★★**本ファイル**（群準同型） |
| — | `PicSheaf (Spec 𝓞_F) ≃* Cl(F)` | ✅ 在庫（`picSheafEquivClassGroupOF`） |

## ★★★★★★★★在庫が全部持っていた

★測定（2026-08-28）で分かったのは、要る部品が**すべて在庫**だったことである:

| 部品 | 場所 |
|---|---|
| `isLocallyTrivial_sheafify`（層化は局所自明性を保つ） | `PicSheafifyTrivial.lean` |
| `InvSheaf.ofLocallyTrivial`（局所自明な層加群は可逆層） | `PicEvalIso.lean` |
| `PicSheaf.mk_eq_mk`（同型なら同じ類） | `PicSheafGroup.lean` |
| `picSheafEquivClassGroupOF`（`Pic(Spec 𝓞_F) ≃* Cl(F)`） | `PicClassGroup.lean` |

★★したがって本ファイルは**それらを継ぐだけ**である。

## ★★★★★★★★群準同型であること —— 層化はテンソル積と可換

`map_mul` に要るのは

    `層化 (A ⊗ B) ≅ tensorModules (層化 A) (層化 B)`

であり、★これも在庫を 2 本継ぐだけであった:

    `sheafifyTensorRight`：`層化 (A ⊗ B) ≅ 層化 ((層化 A).val ⊗ B)`
    `sheafifyTensorLeft` ：`層化 (A' ⊗ B) ≅ 層化 (A' ⊗ (層化 B).val)`

★★`map_one` は要らない——`APicM X` も `PicSheaf X` も**群**なので
`MonoidHom.mk'` が `map_mul` だけから作る。

## ★残っている段（明示）

★★アルキメデス成分（計量から無限素点の実数を読む段）と、
`ADiv/APrc ≅ APicM` の組み立て（逆向きを含む）が残る。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace NumberField
open ABC3.Found.GenEll

/-! ## ★★★★★★層化して可逆層と見る -/

/-- ★★★★★★**計量つき算術直線束を層の水準の可逆層と見る**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★層化して `InvSheaf.ofLocallyTrivial` に渡すだけである
——局所自明性は `isLocallyTrivial_sheafify`（在庫）が保つ。 -/
noncomputable def AInv.toInvSheaf {X : Scheme.{0}} (L : AInv X) : InvSheaf X :=
  InvSheaf.ofLocallyTrivial ((sheafifyFunctor X).obj L.carrier.sheaf)
    (isLocallyTrivial_sheafify X L.carrier.sheaf L.carrier.triv)

/-- ★★★**等長同型なら層化して同じ類になる**。

★等長性は捨てて同型だけを使う——`Pic` は計量を見ないからである。 -/
theorem picSheaf_mk_eq_of_isometric {X : Scheme.{0}} {L M : AInv X}
    (h : Isometric L.carrier M.carrier) :
    PicSheaf.mk (AInv.toInvSheaf L) = PicSheaf.mk (AInv.toInvSheaf M) :=
  h.elim fun φ _ => PicSheaf.mk_eq_mk ((sheafifyFunctor X).mapIso φ)

/-! ## ★★★★★★★★`APicM X → PicSheaf X` -/

/-- ★★★★★★★★**`APic`（前層・計量つき）から `Pic`（層）への写像**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★これが台帳 `arakelov-degF-finite-places` の**段 A3**（写像の水準）である。
★★★群準同型版は下の `APicMToPicSheafHom`。 -/
noncomputable def APicMToPicSheaf {X : Scheme.{0}} : APicM X → PicSheaf X :=
  Quotient.lift (fun L => PicSheaf.mk (AInv.toInvSheaf L))
    (fun _ _ h => picSheaf_mk_eq_of_isometric h)

@[simp] theorem APicMToPicSheaf_mk {X : Scheme.{0}} (L : AInv X) :
    APicMToPicSheaf (APicM.mk L) = PicSheaf.mk (AInv.toInvSheaf L) := rfl

/-! ## ★★★★★★★★★★有限素点側——類群へ -/

/-- ★★★★★★★★★★**`APicM (Spec 𝓞_F) → Cl(F)`**——`deg_F` の橋の有限素点側。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★在庫の `picSheafEquivClassGroupOF`（`Pic(Spec 𝓞_F) ≃* Cl(F)`）を継ぐだけである。
★★群準同型版は下の `APicMToClassGroupHom`。 -/
noncomputable def APicMToClassGroup (F : Type) [Field F] [NumberField F] :
    APicM (Spec (CommRingCat.of (𝓞 F))) → ClassGroup (𝓞 F) :=
  fun L => picSheafEquivClassGroupOF F (APicMToPicSheaf L)

@[simp] theorem APicMToClassGroup_mk (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F)))) :
    APicMToClassGroup F (APicM.mk L)
      = picSheafEquivClassGroupOF F (PicSheaf.mk (AInv.toInvSheaf L)) := rfl

/-! ## ★★★★★★★★層化はテンソル積と可換である -/

/-- ★★★★★★★★**層化はテンソル積と可換である**（可逆層について）。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★在庫の `sheafifyTensorRight` と `sheafifyTensorLeft` を継ぐだけである。 -/
noncomputable def sheafifyTensorIso {X : Scheme.{0}} (A B : X.PresheafOfModules)
    (hA : IsLocallyTrivial X A) (hB : IsLocallyTrivial X B) :
    (sheafifyFunctor X).obj (A ⊗ B)
      ≅ tensorModules ((sheafifyFunctor X).obj A) ((sheafifyFunctor X).obj B) :=
  sheafifyTensorRight X A B hB.isLocallyRankOne
    ≪≫ sheafifyTensorLeft X ((sheafifyFunctor X).obj A).val B
        (isLocallyTrivial_sheafify X A hA).isLocallyRankOne

theorem APicMToPicSheaf_mul_mk {X : Scheme.{0}} (L M : AInv X) :
    APicMToPicSheaf (APicM.mk L * APicM.mk M)
      = APicMToPicSheaf (APicM.mk L) * APicMToPicSheaf (APicM.mk M) := by
  show PicSheaf.mk (AInv.toInvSheaf (L.mul M))
    = PicSheaf.mk ((AInv.toInvSheaf L).mul (AInv.toInvSheaf M))
  exact PicSheaf.mk_eq_mk
    (sheafifyTensorIso L.carrier.sheaf M.carrier.sheaf L.carrier.triv M.carrier.triv)

/-- ★★★★★★★★★★**`APicM X →* PicSheaf X` は群準同型である**——台帳の段 A3。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★`map_one` は要らない——どちらも**群**なので `MonoidHom.mk'` が
`map_mul` だけから作る。 -/
noncomputable def APicMToPicSheafHom {X : Scheme.{0}} : APicM X →* PicSheaf X :=
  MonoidHom.mk' APicMToPicSheaf (by
    refine Quotient.ind fun L => Quotient.ind fun M => ?_
    exact APicMToPicSheaf_mul_mk L M)

/-- ★★★★★★★★★★**`APicM (Spec 𝓞_F) →* Cl(F)`**——`deg_F` の橋の有限素点側、
**群準同型として**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★これで `ADiv(F)/APrc(F) ≅ APic(Spec 𝓞_F)` の**有限素点側が群準同型として存在する**。
★★★残るのはアルキメデス成分と、逆向きの構成である。 -/
noncomputable def APicMToClassGroupHom (F : Type) [Field F] [NumberField F] :
    APicM (Spec (CommRingCat.of (𝓞 F))) →* ClassGroup (𝓞 F) :=
  (picSheafEquivClassGroupOF F).toMonoidHom.comp APicMToPicSheafHom

@[simp] theorem APicMToClassGroupHom_mk (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F)))) :
    APicMToClassGroupHom F (APicM.mk L)
      = picSheafEquivClassGroupOF F (PicSheaf.mk (AInv.toInvSheaf L)) := rfl

/-! ## ★★★★★★★★★★段 A の到達点 —— `Γ(L)` は可逆加群である -/

/-- ★★★★★★★★★★**台帳の段 A（アフィンの橋）そのもの**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

台帳 `arakelov-degF-finite-places` の段 A の claim は

> `Spec R` 上の局所自明な**前層**加群 `L` から**可逆 `R`-加群** `Γ(L,⊤)` を得る

であった。★層化して（`AInv.toInvSheaf`）在庫の `invertible_gammaCarrier` に渡すだけである。

★★これで**段 A が閉じた**。★★★段 B（有限部分 `log #(Γ(L)/𝓞_F·s)`）の前提が揃う
——`Γ(L)` が可逆（＝階数 1 射影）なので、`s ≠ 0` に対して商が有限になる。 -/
theorem invertible_gamma_toInvSheaf (R : CommRingCat.{0}) (L : AInv (Spec R)) :
    Module.Invertible (R : Type)
      (AlgebraicGeometry.moduleSpecΓFunctor.obj (AInv.toInvSheaf L).carrier : Type) :=
  invertible_gammaCarrier R (AInv.toInvSheaf L)

/-! ## ★★★★★段 B への一歩 —— 非零切断の存在 -/

open scoped TensorProduct in
/-- ★★★**可逆加群は自明でない**。

★`Mv ⊗ M ≃ R` なので、`M = 0` なら `R = 0` になる。 -/
theorem invertible_nontrivial (R M : Type) [CommRing R] [AddCommGroup M] [Module R M]
    [Module.Invertible R M] [Nontrivial R] : Nontrivial M := by
  by_contra h
  rw [not_nontrivial_iff_subsingleton] at h
  haveI : Subsingleton (Module.Dual R M ⊗[R] M) := by infer_instance
  haveI : Subsingleton R := (Module.Invertible.linearEquiv R M).symm.toEquiv.subsingleton
  exact false_of_nontrivial_of_subsingleton R

/-- ★★★★★**算術直線束には非零の大域切断がある**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★古典的な `deg_F` の式は **非零の `s`** を取る。
★★段 A が閉じたのでその存在が言える。 -/
theorem exists_ne_zero_gamma (R : CommRingCat.{0}) [Nontrivial (R : Type)]
    (L : AInv (Spec R)) :
    ∃ s : (AlgebraicGeometry.moduleSpecΓFunctor.obj (AInv.toInvSheaf L).carrier : Type),
      s ≠ 0 := by
  haveI := invertible_gamma_toInvSheaf R L
  haveI := invertible_nontrivial (R : Type)
    (AlgebraicGeometry.moduleSpecΓFunctor.obj (AInv.toInvSheaf L).carrier : Type)
  exact exists_ne 0

/-! ## ★★★★段 B の前提 —— 有限生成と射影性 -/

/-- ★★★**`Γ(L)` は有限生成である**。

★mathlib の「An invertible module is finite and projective」が
**インスタンスとして**与えてくれる(実測 2026-08-28)。 -/
theorem finite_gamma (R : CommRingCat.{0}) (L : AInv (Spec R)) :
    Module.Finite (R : Type)
      (AlgebraicGeometry.moduleSpecΓFunctor.obj (AInv.toInvSheaf L).carrier : Type) := by
  haveI := invertible_gamma_toInvSheaf R L
  infer_instance

/-- ★★★**`Γ(L)` は射影である**。★射影＋整域なので**捻れなし**であり、
これが段 B(商の有限性)の出発点になる。 -/
theorem projective_gamma (R : CommRingCat.{0}) (L : AInv (Spec R)) :
    Module.Projective (R : Type)
      (AlgebraicGeometry.moduleSpecΓFunctor.obj (AInv.toInvSheaf L).carrier : Type) := by
  haveI := invertible_gamma_toInvSheaf R L
  infer_instance

/-! ### ★出典の紐付け(`.src`) -/

def AInv.toInvSheaf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(計量つき算術直線束を層の水準の可逆層と見ること)",
    sectionId := "genell-def-1-1-ii" }

def APicMToPicSheaf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(APicM X → PicSheaf X——台帳の段 A3、写像の水準)",
    sectionId := "genell-def-1-1-ii" }

def APicMToClassGroup.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(APicM (Spec 𝓞_F) → Cl(F)——deg_F の橋の有限素点側)",
    sectionId := "genell-def-1-1-ii" }

def APicMToPicSheafHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(APicM X →* PicSheaf X が群準同型であること——台帳の段 A3)",
    sectionId := "genell-def-1-1-ii" }

def APicMToClassGroupHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(APicM (Spec 𝓞_F) →* Cl(F)——deg_F の橋の有限素点側、群準同型)",
    sectionId := "genell-def-1-1-ii" }

def sheafifyTensorIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(層化はテンソル積と可換であること)",
    sectionId := "genell-def-1-1-ii" }

def invertible_gamma_toInvSheaf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(局所自明な前層加群から可逆 R-加群 Γ(L,⊤)——台帳の段 A)",
    sectionId := "genell-def-1-1-ii" }

def finite_gamma.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(Γ(L) が有限生成・射影であること——段 B の前提)",
    sectionId := "genell-def-1-1-ii" }

def finite_gamma.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Module.Invertible から Module.Finite / Module.Projective(インスタンス)"
      (.inMathlib "Module.Invertible") 4,
    .implicitStep
      ("★段 B の残りの道筋(実測 2026-08-28): " ++
       "Γ(L) は有限生成・射影(mathlib のインスタンス)、" ++
       "O_F は有限生成 Z-加群なので Γ(L) も有限生成 Z-加群。" ++
       "★★残るのは「Γ(L)/(O_F·s) が捻れ」であり、" ++
       "それは Γ(L) ⊗ F が 1 次元 F-ベクトル空間で s ≠ 0 がそれを張ることから出る。" ++
       "★★★有限生成＋捻れの Z-加群は有限") 4 ]

def exists_ne_zero_gamma.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(算術直線束に非零の大域切断があること——段 B への一歩)",
    sectionId := "genell-def-1-1-ii" }

def exists_ne_zero_gamma.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "invertible_gamma_toInvSheaf(Γ(L) が可逆 R-加群＝段 A)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.invertible_gamma_toInvSheaf") 4,
    .citation "[mathlib]" "Module.Invertible.linearEquiv(Mv ⊗ M ≃ R)"
      (.inMathlib "Module.Invertible.linearEquiv") 4,
    .implicitStep
      ("★★残っている段 B の核: Γ(L)/(O_F·s) が**有限**であること。" ++
       "★機構は「Γ(L) は有限生成・階数 1 なので商は捻れ、" ++
       "O_F は Dedekind で剰余体が有限だから有限長」である。" ++
       "★★実測(2026-08-28): mathlib に Module.Invertible の API と Ideal.absNorm はあるが、" ++
       "可逆加群から分数イデアルへの橋を自分で繋ぐ必要がある") 4 ]

def APicMToPicSheaf.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "isLocallyTrivial_sheafify(層化は局所自明性を保つ＝段 A1)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.isLocallyTrivial_sheafify") 4,
    .citation "[ABC3]" "InvSheaf.ofLocallyTrivial(局所自明な層加群は可逆層)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.InvSheaf.ofLocallyTrivial") 4,
    .citation "[ABC3]" "picSheafEquivClassGroupOF(Pic(Spec 𝓞_F) ≃* Cl(F))"
      (.inProject "ABC3" "ABC3.Found.Arakelov.picSheafEquivClassGroupOF") 4,
    .implicitStep
      ("★測定(2026-08-28): 要る部品が**すべて在庫**だった——" ++
       "isLocallyTrivial_sheafify / InvSheaf.ofLocallyTrivial / PicSheaf.mk_eq_mk / " ++
       "picSheafEquivClassGroupOF。本ファイルはそれらを継ぐだけである") 4,
    .citation "[ABC3]" "sheafifyTensorRight / sheafifyTensorLeft(層化はテンソル積と可換)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.sheafifyTensorRight") 4,
    .implicitStep
      ("★★map_one は要らない——APicM X も PicSheaf X も**群**なので " ++
       "MonoidHom.mk' が map_mul だけから作る") 4,
    .implicitStep
      ("★★★残っている段: アルキメデス成分(計量から無限素点の実数を読む段)と " ++
       "ADiv/APrc ≅ APicM の組み立て(逆向きを含む)が残る") 4 ]

end ABC3.Found.Arakelov
