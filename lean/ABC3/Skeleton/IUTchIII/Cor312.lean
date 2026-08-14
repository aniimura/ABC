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

## ★★2026-08-14: `sorry` は消えた。**しかしそれは成果ではない。**

一度は「現在の `Interface` の下でこの命題は偽であり `sorry` は埋まらない」と記録した。
その後、原文へ当たり直して**2段階**で条件を輸入した結果、反例は構成できなくなり、
`cor_3_12` は `sorry` 無しで通るようになった:

| 輸入したもの | 出所(物理ページ) | 効いた先 |
|---|---|---|
| `outputLogVolumes` / `outputLogVolumes_eq` / `qLogVol_mem` | p.184、証明 Step (xi-e)(xi-f) | 不等式 |
| `thetaUnion_isCompact` / `logVol_hull_ne_top_of_isCompact` | p.175 + p.31 + p.127 | 有限性 |

★**我々は1つも導出していない。** 下の証明は輸入した仮説を組み合わせただけである。
**我々のモデルの中で Corollary 3.12 は自明になった**——Scholze–Stix の
「trivial, not false」の、我々のモデル上での再現。判定と経緯は
`Check/IUTchIII/Cor312Degenerate.lean`。

**これは原典への判定ではない。** 原文は実際に証明を書いており(p.174–186 の
Step (i)–(xii))、写していないのはその中身の方である。

★**残る不足**: `thetaUnion_isCompact` の原文の根拠は p.175 の
「the **[easily verified]** compactness」だけで、**検証は書かれていない**
(実測: 論文全体で `compact` 33 件、`¹,°𝒰_{j,v_ℚ}` のコンパクト性を確立する箇所は 0 件)。
次に転写すべきは Theorem 3.11, (i), (a) のモノ解析的整構造 `I(…) ⊆ I^ℚ(…)`
(= log-shell。p.31 に「compact, hence of finite log-volume」と明示)と Proposition 3.9。 -/
theorem cor_3_12 (D : PilotObjectData) :
    thetaLogVol D ≠ ⊤ ∧ qLogVol D ≤ thetaLogVol D := by
  refine ⟨D.logVol_hull_ne_top_of_isCompact _ D.thetaUnion_isCompact, ?_⟩
  have hm := D.qLogVol_mem
  rw [D.outputLogVolumes_eq] at hm
  rw [qLogVol, D.qLogVol_eq]
  exact hm

/-- **結論の後半(不等式)だけは証明できる**(2026-08-14)。

原文の証明 Step (xi-e)(xi-f)(物理 p.184)を `Interface` の
`outputLogVolumes_eq` / `qLogVol_mem` として受けた結果、
不等式は原文が言うとおり「then follows formally」——`sorry` 無しで出る。

★**これは成功ではない。** 内容を作ったのではなく、原文の証明の該当段を
仮説として**輸入した**だけである。詳細と判定は
`Check/IUTchIII/Cor312Degenerate.lean`。 -/
theorem cor_3_12_inequality (D : PilotObjectData) : qLogVol D ≤ thetaLogVol D := by
  have hm := D.qLogVol_mem
  rw [D.outputLogVolumes_eq] at hm
  rw [qLogVol, D.qLogVol_eq]
  exact hm

def cor_3_12_inequality.src : Source :=
  { paper := "IUTchIII", pdfPage := 174, item := "Corollary 3.12 (conclusion)",
    sectionId := "cor-3-12-conclusion" }

/-- 原文 Step (xi-e)(xi-f) を `Interface` の仮説として受けているので、
我々の側に残る依存は無い。**「依存が無い」ではなく「輸入済み」**であることに注意。 -/
def cor_3_12_inequality.needs : List ProofObligation :=
  [ .implicitStep
      "Step (xi-e)(xi-f) を Interface の仮説として輸入した(我々は証明していない)" 184,
    .implicitStep
      "原文 p.184 の直前の文にある限定『subject to the condition』と『perhaps only up to some sort of approximation, as a result of various indeterminacies』を写していない。approximation は原文のどこにも量化されていないため、写せば強さを我々が決めることになる。仮説の強化にあたりうる意図的な単純化" 184 ]

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
    .derivation "Step (i)-(xii)——原文は「relatively concrete consequence」と述べるが、実際の導出は p.174-182 にわたる" 174,
    .implicitStep
      "有限性 −|log(Θ)| ∈ R の出所。正則包は、領域が relatively compact なら λ·O 型の有界集合、そうでなければ I^Q 全体(= 対数体積 +∞)と場合分けで定義される(Remark 3.9.5, (i))。したがって有限性には『可能な像の和集合が relatively compact』が要るが、原文はそれを Corollary 3.12 の文脈で確立していない(2026-08-14 実測: relatively compact は定義の箇所にしか現れない)" 127 ]

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
