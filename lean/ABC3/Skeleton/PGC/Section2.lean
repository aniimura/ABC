import ABC3.Skeleton.PGC.Section1
import ABC3.Skeleton.PGC.Section1Cor13
import ABC3.Skeleton.PGC.Section2Defs
import ABC3.Found.PGC.CountableGenerators
import ABC3.Found.PGC.FilteredGroup
import ABC3.Found.PGC.LocalFieldNorm

/-!
# [pGC] §2 — 命題

設定・記号は `ABC3/Skeleton/PGC/Setup.lean`。
`Definition 2.3`(filtered group)は原典が境界外入力を要求しない純粋な定義なので、
`sorry` 無しで `ABC3/Skeleton/PGC/Setup.lean` に直接置いた(G8/G9 は
`theorem`/`lemma` のみを見るので、`structure` はそもそも対象外)。
★2026-09-08 まで置き場は `ABC3/Found/PGC/FilteredGroup.lean` だった——
`Skeleton/PGC/Section3Defs` / `Section4Defs` が `Found` を import せずに立つよう
`Setup.lean` へ降ろした(旧名は `Found/PGC/FilteredGroup.lean` の `export` で残してある)。

構造化: `ResearchPaper/1_Structured/A Version of the Grothendieck Conjecture for p-adic
Local Fields/section-2.html`(PDF 目視確認 2026-09-03、物理 p.4-5)。
-/

namespace ABC3.Skeleton.PGC

open ABC3.Meta ABC3.Interface.PGC ABC3.Found.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ## Proposition 2.1 -/

/-- **[pGC] Proposition 2.1**

原文 (pGC p.4):
> The Γ_K-module K[bar] may be recovered group-theoretically from Γ_K.

## 未解決

証明は独立した proof ブロックを持たず、直前の地の文(p進対数による U_K ⊗ Q_p ≅ K の同型 +
Verlagerung による有限拡大への遷移の両立性)がそのまま論拠になっている。

## 依拠する境界外の結果

- Verlagerung(転送写像)自体は **mathlib に存在する**(`MonoidHom.transfer`、
  `Mathlib/GroupTheory/Transfer.lean:148`)——§1・§2 でこれまで調べた境界外入力の中で
  初めて「公理化不要」と判明した対象。ただし本命題全体の証明には他に p進対数の
  定量評価(境界外、[5])も必要。

## ★★★2026-09-08: 埋まった —— ★**原典の道を 1 つも通らずに閉じた**

★上に挙げた **Verlagerung も p 進対数も、結局 1 度も使っていない**。
実際に通った道は `Found/PGC/` の 3 本:

- `SmoothModelTransport.lean` —— `K̄ ≅ C^∞(Γ_K, K)` から Prop 2.1 への還元(★汎関数形式)
- `NormalBasisFunctional.lean` —— ★`SmoothModelCarrier K ↔ HasCoherentFunctional K`(**同値**)
- `CoherentFunctional.lean` / `CountableGenerators.lean` —— ★**`HasCoherentFunctional` を無条件に構成**

★**Maschke の群環も、Krasner も要らなかった**:
同変な収縮があれば `Ψ := ι∘φ∘ρ + (1 − ι∘ρ)` が同変自己同型という初等的観察だけで足り、
可算性は **`K` 自身の可算稠密部分集合**(`TopologicalSpace.exists_countable_dense`、
`SeparableSpace` はインスタンスで出る)から出た。

★`#print axioms ABC3.Found.PGC.prop_2_1` に `sorryAx` は無い ——
★**依存の連鎖(上記 3 本すべて)が sorry-free であることの証明にもなっている。** -/
theorem prop_2_1 : RecoverableAsAddModule (p := p) (fun K => K.closure) :=
  ABC3.Found.PGC.prop_2_1

def prop_2_1.src : Source :=
  { paper := "pGC", pdfPage := 4, item := "Proposition 2.1", sectionId := "prop-2-1" }

/-- 原文の証明文から抽出した、証明が要求するもの(G6)。★下界。 -/
def prop_2_1.needs : List ProofObligation :=
  [ .folklore
      ("p進対数が U_K(捩れを法として)を K の開部分群に写す(§2 冒頭)。" ++
       "★2026-09-04: 準同型性は一般の K で確立(`Found/PGC/PadicLogMul.lean::padicLog_mul`)。" ++
       "★★同日: 単射性(‖x‖≤1/4 の球上)も一般の K で確立" ++
       "(`Found/PGC/PadicLogInjective.lean::padicLog_injOn`)。" ++
       "★★★同日: 全射性(縮小写像+Banach の不動点定理、exp/log の互逆性を経由しない" ++
       "別ルート)まで確立し、`padicLog K` が半径 1/4 の球からそれ自身への" ++
       "全単射であることを示した(`Found/PGC/PadicLogSurjective.lean::padicLog_bijOn`)。" ++
       "この folklore 入力は sorry 無しで解消") 4,
    .citation "[6] Serre, Local Fields, Chapter VII §8" "Verlagerung(転送写像)"
      (.inMathlib "MonoidHom.transfer") 4,
    .implicitStep "log による U_K⊗Q_p ≅ K と Verlagerung の両立性から Prop 2.1 自体への一段" 4 ]

/-! ## Proposition 2.2 -/

/-- **[pGC] Proposition 2.2**

原文 (pGC p.5):
> Suppose that we are given the following group-theoretic data: the topological group Γ_K,
> together with the indexed filtration Γ_K^v for all v > 0. Then the Γ_K-modules O[scr]_K[bar], and
> K[bar]∧ can be recovered group-theoretically from this group-theoretic data.

## 形式化上の簡略化(逸脱として記録)

原文の O_K̄(K̄ の整数環)・K̄^(K̄ の p進完備化)は、無限次拡大 K̄ 上の付値の構成
(スペクトルノルムの colimit 的延長)を要し、`Found/PGC/LocalFieldNorm.lean` の
スペクトルノルム機構(有限次拡大にのみ適用)をそのままでは使えない。
ここでは両者を**未構築の対象**として抽象化し(`IntKbar`・`CompKbar`、任意の
Γ_K-加群として与える)、Prop 2.1 と同じ `RecoverableAsAddModule` の形で主張する
——具体的な構成は別途 `Found/` の課題として残す。

## 依拠する境界外の結果

- `RamificationFiltration p`(Herbrand の定理、mathlib 不在——
  `Interface/PGC/LocalFieldData.lean`)。
- 上付き↔下付き番号付けの変換([6] Serre, Chapter IV)。
- `Γ_K^0 = I_K`(Corollary 1.3 の系、§1 への直接依存)。 -/
theorem prop_2_2 (_RF : RamificationFiltration p)
    (IntKbar CompKbar : PAdicLocalField p → Type*)
    [∀ K, AddCommGroup (IntKbar K)] [∀ K, DistribMulAction K.absGal (IntKbar K)]
    [∀ K, AddCommGroup (CompKbar K)] [∀ K, DistribMulAction K.absGal (CompKbar K)] :
    RecoverableAsAddModule IntKbar ∧ RecoverableAsAddModule CompKbar := sorry

def prop_2_2.src : Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- 原文の証明文から抽出した、証明が要求するもの(G6)。★下界。 -/
def prop_2_2.needs : List ProofObligation :=
  [ .citation "Interface.PGC.RamificationFiltration" "高次分岐群(上付き番号付け)"
      (.absent "mathlib v4.31.0-rc2 実測: RamificationGroup.lean に上付き番号付けは無い。★2026-09-06 に再測: RamificationGroup.lean にあるのは分解群・惰性群の 4 宣言だけで、re:`upperNumbering|UpperNumbering|ramificationGroup[A-Za-z]*Upper|Herbrand|herbrand|phiDeriv|psiDeriv`→0") 4,
    .citation "[6] Serre, Local Fields, Chapter IV" "上付き↔下付き番号付けの変換"
      (.absent "mathlib v4.31.0-rc2 実測: lowerNumbering/upperNumbering 系の宣言は0件。★2026-09-06 に再測: re:`lowerNumbering|LowerNumbering|upperNumbering|UpperNumbering|herbrand|Herbrand`→0。★★2026-09-06 第 3 の再測でこの記録は半分誤りと分かった —— 下付き分岐群は名前が違うだけで mathlib に在る: `Ideal.inertia` (RingTheory/Ideal/Defs.lean:152)、`AddSubgroup.inertia` (Algebra/Group/Subgroup/Basic.lean:1066)、`AddSubgroup.mem_inertia` は @[simp]、そして `AddSubgroup.subgroupOf_inertia` は rfl(= 原文の「下付きは部分群への移行と両立する」がそのまま在る)。不在なのは上付き番号付けと Herbrand の変換だけ。★教訓: 我々が付けたい名前で引いて 0 件でも、対象が別名で在ることがある") 5,
    .otherPaper "pGC" "Corollary 1.3" 3 ]

end ABC3.Skeleton.PGC
