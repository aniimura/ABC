import ABC3.Interface.IUTchIII.PilotObjects
import Mathlib.Tactic.Linarith

/-!
# [IUTchIII] Corollary 3.12 — 北極星の statement

原典: S. Mochizuki, *Inter-universal Teichmüller Theory III* (May 2020)、物理 p.173–174。
構造化: `ResearchPaper/1_Structured/Inter-universal Teichmuller Theory III/section-3.html`
(PDF 目視確認済み、二重下線は 400 dpi 拡大で確認)。

## この段階で確定していること / していないこと

- 型は付く。**証明は付けていない**(`sorry`)。
- 未構築の基礎は `Interface/IUTchIII/PilotObjects.lean` の `PilotObjectData` で受けている。
  **`Interface` を置いたことは形式化ではない**(PLAN §9)。
- `PilotObjectData` は **Corollary 3.12 が名指しで読み出す分だけ**を受けており、
  「the situation of Theorem 3.11」の全体ではない(その未列挙は構造化の
  `<p class="open">` に記録済み)。

## 非対称性

原文は Θ 側にだけ不定性を掛ける。その非対称性は本ファイルの2つの定義に現れる:

- `thetaLogVol` — 可能な像の**和集合**の**正則包**を測る(`⋃₀` と `hull` が付く)
- `qLogVol`     — **その像**を測る(どちらも付かない)

`logVol` は**同じもの**を使う——原文が同じ多輻的表現の中で両者を比べているため。

## 結論は2つある

原文 (IUTchIII p.174):
> Then it holds that − |log(Θ[ul2])| ∈ R[bb], and
> − |log(Θ[ul2])| ≥ − |log(q[ul2])|

前半(**有限性**)は後半の前提ではなく、**結論の一部**である。
`WithTop ℝ` では `x ≤ ⊤` が常に成り立つので、有限性を落とすと不等式だけが残り、
`thetaLogVol = ⊤` で自明に成立してしまう。ゆえに2つを連言で述べる。
-/

namespace ABC3.Skeleton.IUTchIII

open ABC3.Meta ABC3.Interface.IUTchIII

/-- **`−|log(Θ[ul2])|`** — Θ-パイロット対象の側の対数体積。

原文 (IUTchIII p.173):
> of the holomorphic hull [cf. Remark 3.9.5, (i)]
> of the union of the possible images of a Θ-pilot object [cf. Definition 3.8, (i)],

`⋃₀`(可能な像の和集合)と `hull`(正則包)が付くのが q 側との違い。 -/
noncomputable def thetaLogVol (D : PilotObjectData) : WithTop ℝ :=
  D.logVol (D.hull (⋃₀ D.possibleThetaImages))

def thetaLogVol.src : Source :=
  { paper := "IUTchIII", pdfPage := 173, item := "Corollary 3.12 (Θ-pilot side)",
    sectionId := "cor-3-12-theta" }

/-- **`−|log(q[ul2])|`** — q-パイロット対象の側の対数体積。

原文 (IUTchIII p.174):
> for the procession-normalized mono-analytic log-volume of the image of a
> q-pilot object [cf. Definition 3.8, (i)], relative to the relevant Kummer isomor-
> phisms [cf. Theorem 3.11, (ii)], in the multiradial representation of Theorem

**単一の像**をそのまま測る。`⋃₀` も `hull` も付かない。 -/
noncomputable def qLogVol (D : PilotObjectData) : WithTop ℝ :=
  D.logVol D.qImage

def qLogVol.src : Source :=
  { paper := "IUTchIII", pdfPage := 174, item := "Corollary 3.12 (q-pilot side)",
    sectionId := "cor-3-12-q" }

/-- **[IUTchIII] Corollary 3.12**(北極星)

原文 (IUTchIII p.173):
> Corollary 3.12. (Log-volume Estimates for Θ-Pilot Objects) Suppose
> that we are in the situation of Theorem 3.11.

原文 (IUTchIII p.174):
> Then it holds that − |log(Θ[ul2])| ∈ R[bb], and
> − |log(Θ[ul2])| ≥ − |log(q[ul2])|

## 条件付き形式化

`D : PilotObjectData` を仮説に取る。原文が §3 で所与として使うが、
我々がまだ構成できていないもの(`Interface/IUTchIII/PilotObjects.lean`)。

## 結論の2つの成分

- `thetaLogVol D ≠ ⊤` — 原文の「it holds that `−|log(Θ)| ∈ R`」。
  Θ 側の値域は `ℝ ⋃ {+∞}` なので、これは**主張**である。
- `qLogVol D ≤ thetaLogVol D` — 原文の不等式。

## ★退化について

この statement は、不定性を自明化した `PilotObjectData`(可能な像が
`{qImage}` ただ1つ、`hull = id`)の下で**等号として自明に成り立つ**。
実証は `Check/IUTchIII/Cor312Degenerate.lean`(`sorry` 無し)。
Scholze–Stix の「偽ではなく自明になる」を我々のモデル上で再現したもので、
**原典への判定ではない**——我々の `Interface` が弱いだけかもしれない。

## ★★この `sorry` は埋まらない(2026-08-14 実測)

さらに調べたところ、**現在の `Interface` の下でこの命題は偽である**——
結論の前半(有限性 `thetaLogVol ≠ ⊤`)を破る `PilotObjectData` が構成できる
(`Check.IUTchIII.cor_3_12_refutable_under_current_interface`、`sorry` 無し)。
したがって**この `sorry` は原理的に埋まらない**。埋めようとしないこと。

トリアージ(PLAN §5-2、**既定は①**):

- **① 我々のモデル化の誤り** ← 現時点の分類。`PilotObjectData` は `hull` にも
  `possibleThetaImages` にも条件を課していないので、原文が Theorem 3.11 から
  受け取っているはずの内容を運べていない。
- ② 必要な数学が未構築(mono-theta 環境・log-shell・Frobenioid 等、mathlib に 0 件)。
- ③ 原典側の飛躍 —— **名乗らない**。§5-2 の要件(複数の独立な型設計・falsifier)を
  満たしていない。

**意味すること**: 原文の有限性は Corollary 3.12 の**結論**なので仮説には置けない。
すなわち **Corollary 3.12 の前半は Corollary 3.12 自身の逐語からは出せず、
Theorem 3.11 の中身が要る**。構造化で `<p class="open">` に記録した
「the situation of Theorem 3.11 が未列挙」の、最初の具体的な代償である。 -/
theorem cor_3_12 (D : PilotObjectData) :
    thetaLogVol D ≠ ⊤ ∧ qLogVol D ≤ thetaLogVol D := by
  sorry

def cor_3_12.src : Source :=
  { paper := "IUTchIII", pdfPage := 174, item := "Corollary 3.12 (conclusion)",
    sectionId := "cor-3-12-conclusion" }

/-- 原文の証明文から抽出した、証明が要求するもの(G6)。★**下界**。

原文 p.174–182 の Step (i)–(xii) を数えたのではなく、
Corollary 3.12 の**本文と証明冒頭が名指ししているもの**だけを拾った。 -/
def cor_3_12.needs : List ProofObligation :=
  [ .otherPaper "[IUTchIII]" "Theorem 3.11(multiradial representation・(Ind1)(Ind2)(Ind3))" 153,
    .otherPaper "[IUTchIII]" "Definition 3.8, (i)(Θ-pilot object / q-pilot object)" 173,
    .otherPaper "[IUTchIII]" "Remark 3.9.5, (i)(holomorphic hull)" 173,
    .otherPaper "[IUTchIII]" "Proposition 3.9, (i), (ii)(手続き正規化モノ解析的対数体積)" 173,
    .otherPaper "[IUTchII]" "Corollary 4.10, (i)(ラベル 0 と ⟨F_l⟩ の同定、記号 △)" 174,
    .otherPaper "[IUTchI]" "Definition 3.1, (b)(楕円曲線 E_F の q-パラメータ)" 174,
    .implicitStep "「the situation of Theorem 3.11」が何を含むかは列挙されていない(物理 p.153-159 の7ページ)" 173,
    .derivation "Step (i)-(xii)——原文は「relatively concrete consequence」と述べるが、実際の導出は p.174-182 にわたる" 174 ]

/-- 原文が「i.e.」で言い換える形。**abc へ効くのはこの形**。

原文 (IUTchIII p.174):
> — i.e., C_Θ ≥ −1 for any real number C_Θ ∈ R[bb] such that − |log(Θ[ul2])| ≤ C_Θ · |log(q[ul2])|.

原文が「i.e.」と書いているとおり、これは `cor_3_12` から**実際に直ちに従う**
(下の証明が `sorry` 無しであることが、その「i.e.」に隠れた段が無いことの確認になる)。 -/
theorem cor_3_12_CTheta (D : PilotObjectData) (C : ℝ)
    (h : thetaLogVol D ≤ ((C * D.qAbs : ℝ) : WithTop ℝ)) : -1 ≤ C := by
  have h2 := (cor_3_12 D).2
  rw [qLogVol, D.qLogVol_eq] at h2
  have h3 : ((-D.qAbs : ℝ) : WithTop ℝ) ≤ ((C * D.qAbs : ℝ) : WithTop ℝ) := le_trans h2 h
  have h4 : -D.qAbs ≤ C * D.qAbs := by exact_mod_cast h3
  nlinarith [D.qAbs_pos]

def cor_3_12_CTheta.src : Source :=
  { paper := "IUTchIII", pdfPage := 174, item := "Corollary 3.12 (conclusion)",
    sectionId := "cor-3-12-conclusion" }

/-- 「i.e.」の言い換えは `cor_3_12` から導いたので、外部依存は増えない。 -/
def cor_3_12_CTheta.needs : List ProofObligation := []

end ABC3.Skeleton.IUTchIII
