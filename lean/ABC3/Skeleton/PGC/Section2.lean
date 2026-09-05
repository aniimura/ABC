import ABC3.Skeleton.PGC.Section1
import ABC3.Skeleton.PGC.Section1Cor13
import ABC3.Found.PGC.FilteredGroup
import ABC3.Found.PGC.LocalFieldNorm

/-!
# [pGC] §2 — 命題

設定・記号は `ABC3/Skeleton/PGC/Setup.lean`。
`Definition 2.3`(filtered group)は原典が境界外入力を要求しない純粋な定義なので、
`sorry` 無しで `ABC3/Found/PGC/FilteredGroup.lean` に直接置いた(G8/G9 は
`theorem`/`lemma` のみを見るので、`structure` はそもそも対象外)。

構造化: `ResearchPaper/1_Structured/A Version of the Grothendieck Conjecture for p-adic
Local Fields/section-2.html`(PDF 目視確認 2026-09-03、物理 p.4-5)。
-/

namespace ABC3.Skeleton.PGC

open ABC3.Meta ABC3.Interface.PGC ABC3.Found.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ## Γ_K-加群としての回復可能性(§1 の `RecoverableFromAbsGal` とは別の形)

§1 の3命題(Prop 1.1・1.2・Cor 1.3)はいずれも「単一の値・部分群が回復できる」形
(`AssociatedObject.RecoverableFromAbsGal`)。Proposition 2.1 はこれと異なり、
**加法群としての K̄ の Γ_K-作用込みの構造そのもの**が回復できる、という主張——
「対応する値」ではなく「対応する加法的同型」の存在を述べる。ゆえに新しい形を導入する。 -/

/-- **Γ_K-加群としての回復可能性**。`Obj K` が Γ_K-作用込みの加法群として与えられているとき、
任意の同型 α : Γ_K ≅ Γ_K′ に対し、α と両立する加法的同型 `Obj K ≃+ Obj K′` が存在する
という主張。

## ★★★2026-09-05: 作用のクラスを `SMul` から `DistribMulAction` に直した

原文は「Γ_K-**加群**」と言うが、旧形は作用を **`SMul`**(`one_smul` も
`mul_smul` も `smul_add` も要求しない、公理ゼロのクラス)で受け取っていた。
その形では `prop_2_2` は**偽**である——`Γ_{ℚ_p}` の非可換性
(`Found/PGC/QpNonAbelian.lean`)を使って病的な「作用」
`g • n := if g = g₀ then 0 else n` を作れば反例になる
(`Check/PGC/Prop22Degenerate.lean::prop_2_2_statement_false`、`sorry` 無し)。

`DistribMulAction`(加法的自己同型として作用する)に強めると、この病的な
作用は `MulAction` ですらないので塞がる。`Prop 2.1` が使う `K.closure` への
自然な作用はもちろん `DistribMulAction` を満たす。 -/
def RecoverableAsAddModule (Obj : PAdicLocalField p → Type*)
    [∀ K, AddCommGroup (Obj K)] [∀ K, DistribMulAction K.absGal (Obj K)] : Prop :=
  ∀ {K K' : PAdicLocalField p} (α : ContinuousMulEquiv K.absGal K'.absGal),
    ∃ φ : Obj K ≃+ Obj K', ∀ (g : K.absGal) (x : Obj K),
      φ (g • x) = (α.toMulEquiv g) • (φ x)

/-- §2 冒頭(Proposition 2.1 の準備)で導入する我々自身の定義——原典の項目そのものではない
(bare な `"Proposition 2.1"` と紛れないよう、item 名に注記を付ける)。 -/
def RecoverableAsAddModule.src : Source :=
  { paper := "pGC", pdfPage := 4, item := "Proposition 2.1 (RecoverableAsAddModule)",
    sectionId := "prop-2-1" }

/-- `K.closure`(= K̄)への Γ_K の自然な作用(AlgEquiv としての適用)。

★2026-09-05: `SMul`(公理ゼロ)から `DistribMulAction`(=原文の「Γ_K-加群」)に
強めた。理由は `Check/PGC/Prop22Degenerate.lean` を参照——`SMul` のままだと
`prop_2_2` が**偽**になる。 -/
noncomputable instance closureDistribMulAction (K : PAdicLocalField p) :
    DistribMulAction K.absGal K.closure where
  smul g x := g x
  one_smul _ := rfl
  mul_smul _ _ _ := rfl
  smul_zero g := map_zero g
  smul_add g x y := map_add g x y

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
  定量評価(境界外、[5])も必要。 -/
theorem prop_2_1 : RecoverableAsAddModule (p := p) (fun K => K.closure) := sorry

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
      (.absent "mathlib v4.31.0-rc2 実測: RamificationGroup.lean に上付き番号付けは無い") 4,
    .citation "[6] Serre, Local Fields, Chapter IV" "上付き↔下付き番号付けの変換"
      (.absent "mathlib v4.31.0-rc2 実測: lowerNumbering/upperNumbering 系の宣言は0件") 5,
    .otherPaper "pGC" "Corollary 1.3" 3 ]

end ABC3.Skeleton.PGC
