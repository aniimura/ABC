import ABC3.Skeleton.PGC.Setup
import ABC3.Interface.PGC.LocalFieldData

/-!
# [pGC] §1 — 命題

設定・記号・§1冒頭の暗黙の定義は `ABC3/Skeleton/PGC/Setup.lean`。
未構築の基礎は `ABC3/Interface/PGC/LocalFieldData.lean`。
-/

namespace ABC3.Skeleton.PGC

open ABC3.Meta ABC3.Interface.PGC

variable {p : ℕ} [Fact p.Prime]

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
  mathlib 不在(実測 v4.31.0-rc2)。

## 反証可能性(2026-08-14 監査)

Corollary 3.12 が現在の `Interface` の下で偽だったこと(`Check/IUTchIII/Cor312Degenerate.lean`)を
受けて横断監査した。**この定理は反証できなかった。** 本定理は `Interface` を取らないうえ、
共役による経路は `Check.PGC.cyclotomicCharacter_conj_invariant` で
**閉じていることが証明済み**(円分指標は可換群への準同型なので内部自己同型で動かない)。
探した範囲は `Check/PGC/RefutationAttempts.lean`。 -/
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
      (.absent <|
        "mathlib v4.31.0-rc2 全体を再測定(2026-08-14): " ++
        "grep -rlE 'localTateDuality|local Tate|TateDuality|tateDuality' → 0 件、" ++
        "grep -rlE '\\bTate\\b' → 4 件だが全て別物" ++
        "(EllipticCurve/Weierstrass・RingTheory/Adjoin/Tower・Perfectoid/BDeRham・Perfectoid/FontaineTheta)、" ++
        "grep -rlE 'tateCohomology|TateCohomology|GaloisCohomology' → 0 件、" ++
        "grep -rli 'duality' Mathlib/NumberTheory Mathlib/RepresentationTheory → " ++
        "DirichletCharacter/Orthogonality・MulChar/Duality・Cyclotomic/Galois・Tannaka のみ。" ++
        "公開プロジェクト(ClassFieldTheory・LocalClassFieldTheory)は 2026-08-14 に測定して 0 ヒット" ++
        "(ResearchPaper/lean-ecosystem.json。両リポジトリは本リポジトリに clone されていないので今日は再現していない)") 2,
    .derivation "H^2(K,M) ≅ Z/p^nZ ⟺ M ≅ Z/p^nZ(1) を双対性から導く段" 3,
    .implicitStep "Z_p(1) の同型類が回復できることから、指標 χ が回復できることへの一段" 3 ]


/-! ## Proposition 1.2 -/

/-- Proposition 1.2 が回復されると主張する2つの量を「K に付随する対象」として束ねたもの:
剰余体の元の個数 q と、絶対次数 [K : ℚ_p]。

いずれも**数**なので、移送は恒等——原典の「α によって対応物へ移る」は、数については
「等しい」を意味する。

`q` は `Interface` の `ResidueCardinality` から取る(まだ構築できていないため)。
`[K : ℚ_p]` は mathlib の `Module.finrank` で書ける。 -/
noncomputable def residueCardAndDegreeObject (RD : ResidueCardinality p) :
    AssociatedObject p where
  Obj := fun _ => ℕ × ℕ
  obj := fun K => (RD.card K, Module.finrank ℚ_[p] K.carrier)
  transport := fun _ x => x

def residueCardAndDegreeObject.src : Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.2", sectionId := "prop-1-2" }

/-- **[pGC] Proposition 1.2**

原文 (pGC p.3):
> The number q of elements in the residue field of O[scr]_K, and well as the absolute
> degree [K : Q[bb]_p] of K, can be recovered entirely group-theoretically from Γ_K.

原文の "and well as" は "as well as" の誤植と読めるが、逐語欄では原文どおりに保つ。

## 条件付き形式化

`RD : ResidueCardinality p` を仮説として取る——q を与える構成そのものが未構築だから
(`Interface/PGC/LocalFieldData.lean`)。`RD` が本物になれば、この定理は無条件になる。

## 未解決

原典の論拠は局所類体論の同型 Γ_K^ab ≅ (K^×)^∧ と、そこからの計数
(捩れの prime-to-p 部分が q−1 個、pro-p 商の階数が [K:ℚ_p]+1)。
後者は p進対数を使う——**mathlib にも公開プロジェクトにも無い**(実測)。

## 反証可能性(2026-08-14 監査)

**反証できなかった。** `RD : ResidueCardinality p` が自由なデータであることは
**この定理を反証する役に立たない**——K′ = K の経路は `transport` が恒等なので
自明に閉じ(`Check.PGC.residueCardAndDegree_self_route_closed`)、
残るのは「異なる2体の間の連続同型 α」の構成だけ。しかも
`Check.PGC.refutation_reduces_to_alpha` が示すとおり、その α が1つでも作れれば
`RD` に関係なく次数成分だけで落ちる——すなわち**反証するには本命題自身を偽にするしかない**。
探した範囲は `Check/PGC/RefutationAttempts.lean`。 -/
theorem residueCard_and_degree_recoverable (RD : ResidueCardinality p) :
    (residueCardAndDegreeObject RD).RecoverableFromAbsGal := by
  sorry

def residueCard_and_degree_recoverable.src : Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.2", sectionId := "prop-1-2" }

/-- 原文の証明文から抽出した、証明が要求するもの(G6)。★下界。 -/
def residueCard_and_degree_recoverable.needs : List ProofObligation :=
  [ .citation "[3] Serre, Local Class Field Theory (Cassels-Frohlich, 1967)"
      "相互律による Γ_K^ab ≅ (K^×)^∧"
      (.inProgress "kbuzzard/ClassFieldTheory"
        "抽象的な reciprocityIso はあるが具体形は未完成。LocalCFT/ は2ファイルのみ") 3,
    .folklore "「it is well-known that」0 → U_K → (K^×)^∧ → Ẑ → 0(典拠なし)" 3,
    .citation "p進対数" "U_K の開部分群を(捩れを除いて)K の開部分群と同一視する"
      (.absent <|
        "mathlib v4.31.0-rc2 全体を再測定(2026-08-14): " ++
        "grep -rliE 'padicLog|Padic.log|padic_log|logOneAdd|log_one_add' → " ++
        "6 件すべて実/複素の解析(Analysis/SpecialFunctions/*・Normed/Module/MultipliableUniformlyOn)。" ++
        "Mathlib/NumberTheory/Padics/ 直下 13 ファイルに log 系の定義は無く、" ++
        "'log' の出現は全て Nat.log / Real.log(付値の対数)。" ++
        "'padicExp|expPadic' → 0 件。" ++
        "公開プロジェクト(ClassFieldTheory・LocalClassFieldTheory)は 2026-08-14 に測定して 0 件" ++
        "(ResearchPaper/lean-ecosystem.json。両リポジトリは本リポジトリに clone されていないので今日は再現していない)。" ++
        "★2026-09-04: 級数の収束のみ `Found/PGC/PadicLog.lean` に自前で構築した(sorry 無し、" ++
        "`padicLog`/`logTerm_summable`)——準同型性(log((1+x)(1+y))=log(1+x)+log(1+y))は" ++
        "まだ無いので、この境界外入力は依然として absent のまま残る") 3 ]

end ABC3.Skeleton.PGC
