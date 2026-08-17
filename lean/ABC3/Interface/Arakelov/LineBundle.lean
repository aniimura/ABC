import ABC3.Meta.Claim
import Mathlib.AlgebraicGeometry.IdealSheaf.Functorial
import Mathlib.RingTheory.ClassGroup.Basic
import Mathlib.RingTheory.PicardGroup
import Mathlib.NumberTheory.NumberField.Basic

/-!
# Arakelov 理論のスケルトン(1/3)—— **幾何側: 可逆層と `Pic`**

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★なぜ割るのか

`Interface/GenEll/ArithLineBundle.lean` の `waiting` は
**「スキーム上の直線束(層 B)と `X^arc` 上の hermitian 計量(層 C)」**を
1 本の文字列で待っていた。★**それでは個数が数えられない。**

★★本ファイル群は、その 1 本を**個別に埋められる obligation** に割る。
`tools/check.mjs` の「Interface 実装待ち」の行がそのまま**残作業の件数**になる。

## ★★本ファイルが受ける 3 件(幾何側)

| # | 受けるもの | mathlib の在庫(2026-08-17 実測) |
|---|---|---|
| B1 | `Pic(X)`(可逆層の群) | ★`SheafOfModules.IsLocallyFree` は**ある**が階数が無い。テンソル積は**前層のみ**(`ModuleCat/Presheaf/Monoidal.lean`)、層版は無い |
| B2 | Cartier 因子 → 可逆層 `𝒪_X(D)` | ★`Scheme.IdealSheafData` は**ある**(引き戻し `comap` も。積を保つことは我々が証明済み) |
| B3 | `Pic(Spec 𝓞_F) ≅ ClassGroup 𝓞_F` | ★`ClassGroup` は**ある**。橋が無い |

★★★**B1 が律速である**——層の圏のモノイド構造が要る。
`Mathlib/Algebra/Category/ModuleCat/Presheaf/Monoidal.lean` はあるので、
**層化を通せば届く**位置にはある(2026-08-17 実測)。
-/

namespace ABC3.Interface.Arakelov

open ABC3.Meta AlgebraicGeometry CategoryTheory NumberField

/-! ## ★★B1 —— `Pic(X)` -/

/-- **(B1)** スキーム上の**可逆層の群** `Pic(X)`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★原文の `L` は「`X` 上の直線束」である。★★引き戻しと テンソル積が要る——
高さの加法性(`Proposition 1.4, (i)`)がテンソル積の上で述べられるからである。 -/
structure PicardData where
  /-- `Pic(X)` の台。 -/
  Pic : Scheme.{0} → Type
  /-- テンソル積による群構造。 -/
  group : (X : Scheme.{0}) → CommGroup (Pic X)
  /-- 射に沿った引き戻し `f^* : Pic(Y) → Pic(X)`。 -/
  pullback : {X Y : Scheme.{0}} → (X ⟶ Y) → Pic Y → Pic X
  /-- ★引き戻しは群準同型(テンソル積を保つ)。 -/
  pullback_mul : ∀ {X Y : Scheme.{0}} (f : X ⟶ Y) (L M : Pic Y),
    pullback f (@HMul.hMul _ _ _ (@instHMul _ (group Y).toDivInvMonoid.toMonoid.toMulOneClass.toMul) L M)
      = @HMul.hMul _ _ _ (@instHMul _ (group X).toDivInvMonoid.toMonoid.toMulOneClass.toMul)
          (pullback f L) (pullback f M)
  /-- ★引き戻しは関手的(恒等射)。 -/
  pullback_id : ∀ {X : Scheme.{0}} (L : Pic X), pullback (𝟙 X) L = L
  /-- ★引き戻しは関手的(合成)。 -/
  pullback_comp : ∀ {X Y Z : Scheme.{0}} (f : X ⟶ Y) (g : Y ⟶ Z) (L : Pic Z),
    pullback (f ≫ g) L = pullback f (pullback g L)
  /-- ★★★**アフィンでは mathlib の `Pic R` と一致する**。

  ★★★**これが退化を殺す。**`Pic := PUnit`(自明群)では
  `Pic (Spec ℤ[√-5]) ≃* Pic ℤ[√-5]`(位数 2)が成り立たない。
  ★mathlib の `Mathlib/RingTheory/PicardGroup.lean` に `Pic R` と
  `ClassGroup.equivPic` があるので、これは**書ける条件**である。 -/
  equivPicRing : (R : CommRingCat.{0}) →
    letI := (group (Spec R)); Pic (Spec R) ≃* CommRing.Pic R

/-- Track B は何を作らねばならないか。 -/
def PicardData.waiting : WaitingFor :=
  { what := "(B1) スキーム上の可逆層(局所自由階数 1 の加群層)と、そのテンソル積・引き戻しによる群 Pic(X)"
    trackB := "Found/Arakelov — ★mathlib は `SheafOfModules.IsLocallyFree`(`Algebra/Category/ModuleCat/Sheaf/LocallyFree.lean`)と `Scheme.Modules` の pullback/pushforward 随伴を持つ(2026-08-17 実測)。★★無いのは (a) 階数の概念 (b) **層の圏のテンソル積**——前層版 `ModuleCat/Presheaf/Monoidal.lean` はあるので層化を通せば届く。★これが層 B の律速である" }

/-! ## ★★B2 —— Cartier 因子から可逆層へ -/

/-- **(B2)** 有効 Cartier 因子 `D` から可逆層 `𝒪_X(D)` を作る操作。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★**因子表示と可逆層表示を繋ぐのがこれである。**
我々は `Found/GenEll/` で因子表示の高さを構成した——
本 obligation が埋まれば、それが**可逆層の高さになる**。

★`D` を通らない点という条件が消えるのは、可逆層はいつでも引き戻せるからである。 -/
structure CartierPicData where
  /-- 台となる `Pic`。 -/
  toPicardData : PicardData
  /-- `D ↦ 𝒪_X(D)`。 -/
  ofDivisor : (X : Scheme.{0}) → X.IdealSheafData → toPicardData.Pic X
  /-- ★因子の積は可逆層のテンソル積に移る。 -/
  ofDivisor_mul : ∀ (X : Scheme.{0}) (D E : X.IdealSheafData),
    ofDivisor X (D * E)
      = @HMul.hMul _ _ _
          (@instHMul _ (toPicardData.group X).toDivInvMonoid.toMonoid.toMulOneClass.toMul)
          (ofDivisor X D) (ofDivisor X E)
  /-- ★空因子は単位元。 -/
  ofDivisor_top : ∀ (X : Scheme.{0}),
    ofDivisor X ⊤ = (toPicardData.group X).toDivInvMonoid.toMonoid.toOne.one
  /-- ★★**引き戻しと両立する**——これが高さの底変換不変性を出す段である。 -/
  ofDivisor_pullback : ∀ {X Y : Scheme.{0}} (f : X ⟶ Y) (D : Y.IdealSheafData),
    toPicardData.pullback f (ofDivisor Y D) = ofDivisor X (D.comap f)
  /-- ★★★**`𝒪(D)` が自明になるのは `D` が主因子のときに限る**。

  ★★★**これが `ofDivisor := 1` の退化を殺す。**
  ★`Pic` は (B1) の `equivPicRing` で非自明が強制されているので、
  「全部自明」だと**すべての因子が主因子**になってしまい、
  `Pic (Spec R) ≃* CommRing.Pic R` と矛盾する。

  ★原文の `Definition 1.5, (ii)` が `(−)_red` を Cartier 因子として扱えるのは、
  この対応があるからである。 -/
  IsPrincipalDivisor : (X : Scheme.{0}) → X.IdealSheafData → Prop
  ofDivisor_eq_one_iff : ∀ (X : Scheme.{0}) (D : X.IdealSheafData),
    ofDivisor X D
        = (toPicardData.group X).toDivInvMonoid.toMonoid.toOne.one
      ↔ IsPrincipalDivisor X D
  /-- ★★主因子は積で閉じている(主因子のなす部分モノイド)。 -/
  isPrincipalDivisor_mul : ∀ (X : Scheme.{0}) (D E : X.IdealSheafData),
    IsPrincipalDivisor X D → IsPrincipalDivisor X E → IsPrincipalDivisor X (D * E)
  /-- ★★★**`Spec R` では主因子は単項イデアルに対応する**——意味を固定する。

  ★これが無いと `IsPrincipalDivisor := True` で逃げられる。 -/
  isPrincipalDivisor_affine : ∀ (R : CommRingCat.{0}) (D : (Spec R).IdealSheafData),
    IsPrincipalDivisor (Spec R) D ↔
      ((Scheme.IdealSheafData.equivOfIsAffine D).IsPrincipal)

def CartierPicData.waiting : WaitingFor :=
  { what := "(B2) 有効 Cartier 因子 D から可逆層 𝒪_X(D) を作る操作と、それが積・引き戻しと両立すること"
    trackB := "Found/Arakelov — ★`Scheme.IdealSheafData` と `comap` は mathlib にあり、**引き戻しが積を保つこと**(`comap_mul`)は我々が `Found/GenEll/ComapMul.lean` で証明済み(2026-08-17)。★★残るのは (B1) の可逆層そのものであり、本 obligation は (B1) に従属する" }

/-! ## ★★B3 —— `Spec 𝓞_F` 上での `Pic` -/

/-- **(B3)** `Pic(Spec 𝓞_F) ≅ ClassGroup 𝓞_F`。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

★★★**高さの定義が `Spec 𝓞_F` へ引き戻した先で完結するのは、これがあるからである。**
★原文は `ADiv(F)/APrc(F) ⥲ APic(Spec 𝓞_F)` の形で使う——
その有限素点側が本 obligation である。 -/
structure PicSpecData where
  /-- 台となる `Pic`。 -/
  toPicardData : PicardData
  /-- `Pic(Spec 𝓞_F) ≃ ClassGroup 𝓞_F`。 -/
  equivClassGroup : (F : Type) → [Field F] → [NumberField F] →
    toPicardData.Pic (Spec (CommRingCat.of (𝓞 F))) ≃ ClassGroup (𝓞 F)
  /-- ★★★**その同型は (B1) の `equivPicRing` と整合する**。

  ★★★**これが自由な posit を殺す。**mathlib には
  `ClassGroup.equivPic : ClassGroup R ≃* CommRing.Pic R` があるので、
  本条件は `equivClassGroup` を**完全に決めてしまう**——
  すなわち **(B3) は独立の難所ではなく、(B1) の系である**。 -/
  equivClassGroup_compat : ∀ (F : Type) [Field F] [NumberField F]
    (L : toPicardData.Pic (Spec (CommRingCat.of (𝓞 F)))),
    ClassGroup.equivPic (𝓞 F) (equivClassGroup F L)
      = toPicardData.equivPicRing (CommRingCat.of (𝓞 F)) L

def PicSpecData.waiting : WaitingFor :=
  { what := "(B3) Pic(Spec 𝓞_F) ≅ ClassGroup 𝓞_F —— 数体の整数環上の可逆層が類群で書けること"
    trackB := "Found/Arakelov — ★`ClassGroup` は mathlib にある(`RingTheory/ClassGroup/Basic.lean`)。★★無いのは Pic 側(B1)と、その間の橋である。★Dedekind 環では可逆イデアル ≅ 可逆層なので、(B1) が入れば機械的である" }

/-! ## ★出典の紐付け(`.src`) -/

def PicardData.src : Source :=
  { paper := "GenEll", pdfPage := 3, item := "Definition 1.1, (i)(層 B——可逆層と Pic のみ)",
    sectionId := "genell-def-1-1-i" }

def CartierPicData.src : Source :=
  { paper := "GenEll", pdfPage := 3, item := "Definition 1.1, (i)(層 B——Cartier 因子から可逆層へ)",
    sectionId := "genell-def-1-1-i" }

def PicSpecData.src : Source :=
  { paper := "GenEll", pdfPage := 4, item := "Definition 1.2, (i)(層 B——Spec 𝓞_F 上の Pic)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Interface.Arakelov
