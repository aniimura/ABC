import ABC3.Meta.Claim
import Mathlib.AlgebraicGeometry.IdealSheaf.Functorial
import Mathlib.AlgebraicGeometry.Morphisms.Flat
import Mathlib.RingTheory.ClassGroup.Basic
import Mathlib.RingTheory.PicardGroup
import Mathlib.NumberTheory.NumberField.Basic
import Mathlib.AlgebraicGeometry.Modules.Sheaf
import Mathlib.Algebra.Category.ModuleCat.Presheaf.Monoidal
import Mathlib.Algebra.Category.ModuleCat.Presheaf.Sheafification

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

/-! ## ★★★退化封じの道具 —— 層加群のテンソル積(mathlib だけで書く)

★★★**2026-08-17 に穴が見つかった。** 当初の `PicardData` は

    Pic X := CommRing.Pic Γ(X, ⊤)

を**通してしまう**。これは `equivPicRing` も関手性もすべて満たすが、
★**非アフィンな `X` では数学的に誤り**である(例: `Pic(ℙ¹) = ℤ` だが
`Γ(ℙ¹, 𝒪) = k` なので `CommRing.Pic k = 0`)。

★★構造がアフィンでしか `Pic` を縛っていなかったのが原因である。
これは C1 で見つけた `evalAffine` の穴と**同じ型の見落とし**
——「一部で固定しても、残りが自由なら誤った witness が通る」。

★塞ぎ方: `Pic X` の元に**下にある可逆層**を持たせ、
「`Pic X` は可逆層の同型類そのものである」ことを課す。
★`Interface/` は `Found/` を import できないので、テンソル積を
mathlib の前層テンソル(`PresheafOfModules.Monoidal.tensorObj`)+ 層化で書き下す。 -/

/-- ★層加群のテンソル積(前層でテンソルしてから層化する)。 -/
noncomputable def modTensor (X : Scheme.{0}) (F G : X.Modules) : X.Modules :=
  (PresheafOfModules.sheafification (R := X.ringCatSheaf) (𝟙 X.ringCatSheaf.obj)).obj
    (PresheafOfModules.Monoidal.tensorObj (R := X.presheaf) F.val G.val)

/-- ★★**可逆層(直線束)** —— 「テンソル積の逆を持つ」かつ「局所的に階数 1 自由」。

★★★**2 条を課すのは、どちらも直線束の性質であり、どちらも使うからである**:

| 条 | 何に効くか |
|---|---|
| 逆の存在 | `Pic X` の**逆元** |
| 局所自由性 | 結合律の**局所論法**(`Found/Arakelov/PicLocalBasis.lean`) |

★数学的には 2 条は同値である(片方から他方は Nakayama 型の議論で出る)。
★★したがって**内容は変わらない**——Lean で両方を使えるようにしただけである。

★局所自由性は**層論の標準の定義**で述べる——「各開集合が、
`F|_V ≅ 𝒪_V`(**層として**)となる開集合で覆われる」。

★★★**2026-08-18 に強めた。**それまでは「切断ごとに `𝒪(V) ≃ₗ F(V)`」と
書いてあったが、これは標準の定義より**弱い**(制限射との両立を要求しない)。
`Found` 側の `IsLocallyTrivial` はもともと層論の定義で書いてあり、
弱い版のままでは `sheafOf_surjective` が**過大な要求**になっていた。 -/
def IsInvertibleSheaf {X : Scheme.{0}} (F : X.Modules) : Prop :=
  (∃ G : X.Modules, Nonempty (modTensor X F G ≅ SheafOfModules.unit X.ringCatSheaf)) ∧
  (∀ U : X.Opens, ∃ S : Sieve U, S ∈ (Opens.grothendieckTopology X) U ∧
    ∀ ⦃V : X.Opens⦄ (i : V ⟶ U), S i →
      Nonempty ((PresheafOfModules.pushforward₀OfCommRingCat
            (Over.forget V) X.presheaf).obj F.val
          ≅ MonoidalCategory.tensorUnit (C := PresheafOfModules.{0}
            (((Over.forget V).op ⋙ X.presheaf) ⋙ forget₂ CommRingCat.{0} RingCat.{0}))))

/-- **(B1)** スキーム上の**可逆層の群** `Pic(X)`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★原文の `L` は「`X` 上の直線束」である。★★引き戻しと テンソル積が要る——
高さの加法性(`Proposition 1.4, (i)`)がテンソル積の上で述べられるからである。 -/
structure PicardData where
  /-- `Pic(X)` の台。 -/
  Pic : Scheme.{0} → Type 1
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
  /-- ★★★★**`Pic X` の元は可逆層である**——下にある層を持たせる。

  ★★★**これが「非アフィンで自由」の穴を塞ぐ。**
  これが無いと `Pic X := CommRing.Pic Γ(X, ⊤)` が通ってしまい、
  `Pic(ℙ¹) = ℤ` を `0` と主張する witness が「達成」と数えられる。 -/
  sheafOf : (X : Scheme.{0}) → Pic X → X.Modules
  /-- ★下にある層は可逆である。 -/
  sheafOf_invertible : ∀ (X : Scheme.{0}) (L : Pic X), IsInvertibleSheaf (sheafOf X L)
  /-- ★単位元の下にあるのは構造層。 -/
  sheafOf_one : ∀ (X : Scheme.{0}),
    Nonempty (sheafOf X (group X).toDivInvMonoid.toMonoid.toOne.one
      ≅ SheafOfModules.unit X.ringCatSheaf)
  /-- ★★積はテンソル積に移る。 -/
  sheafOf_mul : ∀ (X : Scheme.{0}) (L M : Pic X),
    Nonempty (sheafOf X (@HMul.hMul _ _ _
        (@instHMul _ (group X).toDivInvMonoid.toMonoid.toMulOneClass.toMul) L M)
      ≅ modTensor X (sheafOf X L) (sheafOf X M))
  /-- ★★★★★**引き戻しは層の引き戻しである**。

  ★★★**2026-08-18 に 2 つ目の穴が見つかった。**
  上の `pullback` / `pullback_mul` / `pullback_id` / `pullback_comp` は
  **`sheafOf` と一切結ばれていない**——`pullback` が
  「層の引き戻し」であることをどの条件も要求していなかった。

  ★したがって `pullback f := 1`(自明準同型、ただし恒等射のときだけ `id`)のような
  **幾何と無関係な witness** が通ってしまう。
  ★★これは C1 の `evalAffine`、B1 の `sheafOf` と**同じ型の見落とし**である
  ——「一部で固定しても、残りが自由なら誤った witness が通る」。

  ★★★塞ぎ方: 引き戻しの下にある層が、mathlib の
  `Scheme.Modules.pullback`(層加群の引き戻し)と一致することを課す。 -/
  sheafOf_pullback : ∀ {X Y : Scheme.{0}} (f : X ⟶ Y) (L : Pic Y),
    Nonempty (sheafOf X (pullback f L)
      ≅ (Scheme.Modules.pullback f).obj (sheafOf Y L))
  /-- ★★★**同型な層は同じ元を与える**(`Pic` は同型類の集合)。 -/
  sheafOf_injective : ∀ (X : Scheme.{0}) (L M : Pic X),
    Nonempty (sheafOf X L ≅ sheafOf X M) → L = M
  /-- ★★★**可逆層はすべて現れる**(`Pic` に足りない元が無い)。

  ★★これと `sheafOf_injective` を合わせて「`Pic X` = 可逆層の同型類」が確定する。 -/
  sheafOf_surjective : ∀ (X : Scheme.{0}) (F : X.Modules), IsInvertibleSheaf F →
    ∃ L : Pic X, Nonempty (sheafOf X L ≅ F)

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
  /-- ★★★**`D` が Cartier(可逆イデアル層)であること**。

  ★★★★**2026-08-18 の修正。**それ以前は `ofDivisor` を `IdealSheafData` 全体に
  課していたが、その形は**充足不可能**であった:

  `R = ℚ[x,y]` は UFD なので `CommRing.Pic R = 1`、B1 の `equivPicRing` から
  `Pic(Spec R) = 1`。すると旧 `ofDivisor_eq_one_iff`(無条件)が
  `IsPrincipalDivisor (Spec R) D` を**すべての** `D` に強制し、
  `isPrincipalDivisor_affine` から `R` のすべてのイデアルが単項になる。
  `(x,y)` は単項でないので矛盾する。

  ★原因: `𝒪_X(D)` は **Cartier 因子にしか定義できない**。
  ★★以降の欄はすべて `IsCartierDivisor` の仮定の下で述べる。 -/
  IsCartierDivisor : (X : Scheme.{0}) → X.IdealSheafData → Prop
  /-- ★空因子は Cartier。 -/
  isCartierDivisor_top : ∀ (X : Scheme.{0}), IsCartierDivisor X ⊤
  /-- ★Cartier は積で閉じる。 -/
  isCartierDivisor_mul : ∀ (X : Scheme.{0}) (D E : X.IdealSheafData),
    IsCartierDivisor X D → IsCartierDivisor X E → IsCartierDivisor X (D * E)
  /-- ★Cartier は**平坦**な引き戻しで保たれる。

  ★★★★**2026-08-19 の修正。**それ以前は平坦性を課していなかったが、その形は
  **偽**である:

  `Y = Spec k[x]`、`D = (x)`(`x` は非零因子なので `(x) ≅ k[x]`、可逆)、
  `f : Spec k ⟶ Spec k[x]` を原点とすると、`D.comap f`(= 概型論的逆像の定義イデアル)は
  `Spec k ×_{Spec k[x]} Spec k = Spec k` の恒等射の核、すなわち `⊥` である。
  `Module.Invertible k (⊥ : Ideal k)` は**偽**(階数 0)。

  ★★★古典的にも「Cartier 因子の引き戻しが Cartier」は `f(X) ⊄ Supp D` か
  `f` が平坦のときに限る。可逆**層** `𝒪(D)` の引き戻しは常に可逆だが、
  イデアル層の引き戻しは可逆とは限らない。

  ★下流(高さの底変換不変性)が使うのは `Spec 𝓞_L ⟶ Spec 𝓞_K` で、これは平坦なので
  平坦性を課しても影響しない。 -/
  isCartierDivisor_comap : ∀ {X Y : Scheme.{0}} (f : X ⟶ Y) [Flat f] (D : Y.IdealSheafData),
    IsCartierDivisor Y D → IsCartierDivisor X (D.comap f)
  /-- ★★★**アフィンでは可逆イデアルに対応する**——`IsCartierDivisor` の意味を固定する。

  ★これが無いと `IsCartierDivisor := fun _ _ => False` で全欄が空虚に成立してしまう。 -/
  isCartierDivisor_affine : ∀ (R : CommRingCat.{0}) (D : (Spec R).IdealSheafData),
    IsCartierDivisor (Spec R) D ↔
      Module.Invertible (Γ(Spec R, ⊤) : Type)
        ((Scheme.IdealSheafData.equivOfIsAffine D : Ideal (Γ(Spec R, ⊤) : Type)) : Type)
  /-- `D ↦ 𝒪_X(D)`。 -/
  ofDivisor : (X : Scheme.{0}) → X.IdealSheafData → toPicardData.Pic X
  /-- ★因子の積は可逆層のテンソル積に移る(Cartier のとき)。 -/
  ofDivisor_mul : ∀ (X : Scheme.{0}) (D E : X.IdealSheafData),
    IsCartierDivisor X D → IsCartierDivisor X E →
    ofDivisor X (D * E)
      = @HMul.hMul _ _ _
          (@instHMul _ (toPicardData.group X).toDivInvMonoid.toMonoid.toMulOneClass.toMul)
          (ofDivisor X D) (ofDivisor X E)
  /-- ★空因子は単位元。 -/
  ofDivisor_top : ∀ (X : Scheme.{0}),
    ofDivisor X ⊤ = (toPicardData.group X).toDivInvMonoid.toMonoid.toOne.one
  /-- ★★**平坦な引き戻しと両立する**——これが高さの底変換不変性を出す段である。

  ★平坦性は `isCartierDivisor_comap` と同じ理由(2026-08-19 の修正)。 -/
  ofDivisor_pullback : ∀ {X Y : Scheme.{0}} (f : X ⟶ Y) [Flat f] (D : Y.IdealSheafData),
    IsCartierDivisor Y D →
    toPicardData.pullback f (ofDivisor Y D) = ofDivisor X (D.comap f)
  /-- ★★★**`𝒪(D)` が自明になるのは `D` が主因子のときに限る**(Cartier のとき)。

  ★★★**これが `ofDivisor := 1` の退化を殺す。**
  ★`Pic` は (B1) の `equivPicRing` で非自明が強制されているので、
  「全部自明」だと**すべての可逆イデアルが単項**になってしまい、
  `ℤ[√-5]` の `(2, 1+√-5)`(可逆・非単項)と矛盾する。

  ★原文の `Definition 1.5, (ii)` が `(−)_red` を Cartier 因子として扱えるのは、
  この対応があるからである。 -/
  IsPrincipalDivisor : (X : Scheme.{0}) → X.IdealSheafData → Prop
  ofDivisor_eq_one_iff : ∀ (X : Scheme.{0}) (D : X.IdealSheafData),
    IsCartierDivisor X D →
    (ofDivisor X D
        = (toPicardData.group X).toDivInvMonoid.toMonoid.toOne.one
      ↔ IsPrincipalDivisor X D)
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
    trackB := "Found/Arakelov — ★**(B1) は 2026-08-18 に達成済み**(`Found/Arakelov/PicWitness.lean`)。★★mathlib 在庫を実測(2026-08-18): `Scheme.IdealSheafData` / `Mul` / `comap` / `equivOfIsAffine` は**有る**が、**Cartier 因子は 0 件**(`grep CartierDivisor|EffectiveCartier` → 該当なし)。★★★したがって「可逆イデアル層」の定義から自前で作る。`comap_mul` は `Found/GenEll/ComapMul.lean` に証明済み(2026-08-17)" }

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
