import ABC3.Meta.Claim
import Mathlib.FieldTheory.KrullTopology
import Mathlib.NumberTheory.Cyclotomic.CyclotomicCharacter
import Mathlib.NumberTheory.Padics.PadicNumbers
import Mathlib.FieldTheory.IsAlgClosed.AlgebraicClosure
import Mathlib.Topology.Algebra.ContinuousMonoidHom

/-!
# [pGC] §1 — The Cyclotomic Character and Inertia

原典: S. Mochizuki, *A Version of the Grothendieck Conjecture for p-adic Local Fields* (1997)。
構造化: `ResearchPaper/1_Structured/A Version of the Grothendieck Conjecture for p-adic Local Fields/section-1.html`
(PDF 目視確認済み: 物理 p.2・p.3・p.4冒頭)。

## この段階で確定していること / していないこと

- 型は付く。**証明は付けていない**(`sorry`)。
- §1 冒頭の「group-theoretically に回復できる」は、原典が**見出し語を持たない地の文**で
  与えている定義。ここで明示的に切り出した(`RecoverableFromAbsGal`)。
- ★**原典が暗黙にしていて、我々が明示にしたもの**: 対象を α に沿って移す規則(`transport`)。
  原典は「the object associated to K is necessarily taken by α to the corresponding object
  associated to K′」と述べるだけで、移送の規則自体は文脈から自明として与えていない。
  これをデータとして持たせたのは**我々の設計判断**であり、原典の転写ではない。
  移送の選び方によって主張の強さが変わるため、下記の識別力検査を併置する。
-/

namespace ABC3.Skeleton.PGC

open ABC3.Meta

variable (p : ℕ) [Fact p.Prime]

/-! ## 設定 -/

/-- 本論文でいう **p進局所体** — `ℚ_[p]` の有限次拡大。

原典 [pGC] 物理 p.2:
「Let p be a prime number. Let K be a p-adic local field. By this, we shall mean in this
paper that K is a finite extension of **Q**_p.」 -/
structure PAdicLocalField where
  carrier : Type
  [isField : Field carrier]
  [isAlgebra : Algebra ℚ_[p] carrier]
  [isFinite : FiniteDimensional ℚ_[p] carrier]

def PAdicLocalField.src : Source :=
  { paper := "pGC", pdfPage := 2, item := "Section 1 (opening)", sectionId := "setup-1-intro" }

attribute [instance] PAdicLocalField.isField PAdicLocalField.isAlgebra PAdicLocalField.isFinite

namespace PAdicLocalField

variable {p}

/-- 原典の `K̄` — K の代数閉包(原典は「Fix an algebraic closure K̄ of K」)。 -/
abbrev closure (K : PAdicLocalField p) : Type := AlgebraicClosure K.carrier

def closure.src : Source :=
  { paper := "pGC", pdfPage := 2, item := "Section 1 (opening)", sectionId := "setup-1-intro" }

/-- 原典の `Γ_K ≝ Gal(K̄/K)` — 絶対 Galois 群。

位相は mathlib の `krullTopology` が与える(`Mathlib/FieldTheory/KrullTopology.lean`)。
原典は「isomorphism of topological groups」と言うので、位相込みで扱う必要がある。 -/
abbrev absGal (K : PAdicLocalField p) : Type := K.closure ≃ₐ[K.carrier] K.closure

def absGal.src : Source :=
  { paper := "pGC", pdfPage := 2, item := "Section 1 (opening)", sectionId := "setup-1-intro" }

end PAdicLocalField

/-! ## §1 冒頭の暗黙の定義 -/

/-- **K に付随する対象**の族と、Γ_K の同型に沿った移送の規則。

原典は「an object associated to K」という言い方をするだけで、対象の種類を限定していない
——実際 §1 では3種類が現れる: 指標 χ(Γ_K 上の関数)・数 q や [K:ℚ_p](Γ_K に依存しない値)・
惰性群 I_K(Γ_K の部分群)。ゆえに対象の**型そのもの**が K ごとに変わりうる形にしてある。

★`transport` は原典が明示していない(上記モジュール docstring 参照)。 -/
structure AssociatedObject where
  /-- K に付随する対象の型 -/
  Obj : PAdicLocalField p → Type
  /-- K に付随する対象そのもの -/
  obj : (K : PAdicLocalField p) → Obj K
  /-- 位相群の同型 α : Γ_K ≅ Γ_K′ に沿った移送 -/
  transport : {K K' : PAdicLocalField p} →
    ContinuousMulEquiv K.absGal K'.absGal → Obj K → Obj K'

def AssociatedObject.src : Source :=
  { paper := "pGC", pdfPage := 2, item := "Section 1 (opening)", sectionId := "setup-1-intro" }

variable {p}

/-- **「Γ_K から group-theoretically に回復できる」** — [pGC] §1 冒頭の暗黙の定義。

原文 (pGC p.2):
> when we say that an object associated to K can be recovered "group-theoretically" from
> Γ_K, we mean that given another local p-adic field K[prime], together with an isomorphism
> of topological groups α : Γ_K ≅ Γ[prime]_K, the object associated to K is necessarily
> taken by α to the corresponding object associated to K[prime].

注意すべき形: これは**単称の構成手続きの存在**ではなく、**すべての同型に対する両立性**である。
「Γ_K から χ を作るアルゴリズムがある」ではなく「どの同型で写しても対応が保たれる」。 -/
def AssociatedObject.RecoverableFromAbsGal (A : AssociatedObject p) : Prop :=
  ∀ {K K' : PAdicLocalField p} (α : ContinuousMulEquiv K.absGal K'.absGal),
    A.transport α (A.obj K) = A.obj K'

def AssociatedObject.RecoverableFromAbsGal.src : Source :=
  { paper := "pGC", pdfPage := 2, item := "Section 1 (opening)", sectionId := "setup-1-intro" }

/-!  ## 識別力の検査

`RecoverableFromAbsGal` が常に真でも常に偽でもないことの検査は、主語が違う
(原典の主張ではなく**我々のモデル**についての事実)ため `ABC3/Check/PGC/Section1Discriminating.lean`
に置いた。`Skeleton/` から `Check/` を import してはならない。
-/

/-! ## Proposition 1.1 -/

/-- 円分指標 `χ_K : Γ_K → ℤ_[p]ˣ` を「K に付随する対象」として束ねたもの。

`χ_K` 自体は mathlib の `cyclotomicCharacter` から得る——**これは K(と K̄)を使って定義される**。
Proposition 1.1 が主張するのは、こうして作った対象が実は Γ_K だけで決まる、ということ。
移送は「α で引き戻す」= `f ∘ α.symm`。 -/
noncomputable def cyclotomicCharacterObject : AssociatedObject p where
  Obj := fun K => K.absGal → ℤ_[p]ˣ
  obj := fun K g => cyclotomicCharacter K.closure p g.toRingEquiv
  transport := fun α f g' => f (α.toMulEquiv.symm g')

def cyclotomicCharacterObject.src : Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **[pGC] Proposition 1.1**

原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.

## 未解決(1_Structured の `<p class="open">` と対応)

原典は Proposition 1.1 に独立した証明を付けておらず、直前の地の文(双対性の復習)が論拠になる。
しかしその地の文が確立するのは「Γ_K-加群 **Z**_p(1) の**同型類**が回復できる」ことであり、
本命題が主張するのは「**指標** χ が回復できる」こと。**この間の一段は原典にない。**
証明を付ける段階で、
  ① 我々のモデル化の誤り(§1冒頭の定義から直ちに従う)か
  ② 必要な数学が未構築(局所 Tate 双対性——mathlib 不在)か
  ③ 原典側の飛躍か
を切り分けること。既定は①。

## 依拠する境界外の結果(Track B のキュー候補)

- [2] Bloch–Kato, Grothendieck Festschrift I (1990), Proposition 3.8 — 局所 Tate 双対性。
  mathlib 不在(実測 v4.31.0-rc2)。 -/
theorem cyclotomicCharacter_recoverable :
    (cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal := by
  sorry

def cyclotomicCharacter_recoverable.src : Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- 原文の証明文から抽出した、証明が要求するもの(G6)。

★これは **下界** である——原文が挙げていない依存は写らない。 -/
def cyclotomicCharacter_recoverable.needs : List ProofObligation :=
  [ .citation "[2] Bloch-Kato, Grothendieck Festschrift I (1990)"
      "Proposition 3.8(局所 Tate 双対性)"
      .absent 2,
    .derivation "H^2(K,M) ≅ Z/p^nZ ⟺ M ≅ Z/p^nZ(1) を双対性から導く段" 3,
    .implicitStep "Z_p(1) の同型類が回復できることから、指標 χ が回復できることへの一段" 3 ]

end ABC3.Skeleton.PGC
