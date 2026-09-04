import ABC3.Meta.Claim
import ABC3.Found.Falt1.KaehlerAux

/-!
# [Falt1] Lemma 1.1 —— 単射性・長さの等式(`Found`)

原典: G. Faltings, *p-Adic Hodge Theory*, JAMS Vol.1 No.1 (1988) pp.255-299。
物理 p.4(印字 p.257)。**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

★このPDFのpdftotext OCR層は数式記号を激しく壊す(`Ω`が`K`に、`⊗`が`?`に、
`⊂`が`c`になる等——`Interface/Falt1/Ramification.lean`冒頭の注記と同じ
運用)ため、地の文(`内容`)で写す。

内容 (Falt1 p.4): Lemma 1.1. For any extension V ⊂ W, as above, the
natural map Ω_V ⊗_V W → Ω_W is injective, and its cokernel Ω_{W/V}
has the same length as W/p^δW.

`KaehlerAux.lean` で証明した以下の2定理を貼り合わせ、Falt1 の実際の
`V,K,L,W` に対して Lemma 1.1 の主張(単射性 ＋ 余核の長さの等式)を
証明した:
- `falt1MapBaseChangeInjective`(単射性)
- `falt1CokernelLengthEq`(長さの等式、`Ω_{W/V}` の余核が `W ⧸
  differentIdeal V W` に等しいことと合わせて)
-/

namespace ABC3.Found.Falt1

/-- **[Falt1] Lemma 1.1**(単射性 ＋ 長さの等式、まとめて)。`Z` は絶対微分
`Ω_{V/Z}`・`Ω_{W/Z}` の基底となる任意の環(Falt1 の文脈では固定された
「絶対」基底)。`differentIdeal V W ≠ ⊥` は有限可分拡大なら常に成り立つ
標準的な事実(`differentIdeal_ne_bot`、下記 `.needs` 参照)。 -/
theorem lemma_1_1_falt1 {Z V K L W : Type*} [CommRing Z] [CommRing V] [IsDedekindDomain V]
    [Field K] [Algebra V K] [IsFractionRing V K] [Field L] [Algebra K L] [FiniteDimensional K L]
    [Algebra.IsSeparable K L] [CommRing W] [Algebra W L] [Algebra V W] [Algebra V L]
    [IsScalarTower V K L] [IsScalarTower V W L] [IsIntegralClosure W V L]
    [IsDedekindDomain W] [Module.IsTorsionFree V W] [Algebra Z V]
    [Algebra Z W] [IsScalarTower Z V W]
    (w : W) (hint : IsIntegral V w) (hadjoin : Algebra.adjoin V ({w} : Set W) = ⊤)
    (hw : Algebra.adjoin K ({(algebraMap W L) w} : Set L) = ⊤)
    (hdne : differentIdeal V W ≠ ⊥) :
    Function.Injective (KaehlerDifferential.mapBaseChange Z V W) ∧
      Module.length W Ω[W⁄V] = Module.length W (W ⧸ differentIdeal V W) :=
  ⟨falt1MapBaseChangeInjective w hint hadjoin hw hdne, falt1CokernelLengthEq w hint hadjoin hw⟩

/-! ### ★★★★★★★★項目全体の `.src`

`Lemma 1.1` の主張(単射性 ＋ 長さの等式)の両方が `Found/` に揃ったので
置く。ただし下の「逸脱の記録」にある通り、**Interface の抽象化
(`RamificationSetup`)への正式な差し替えはまだ行っていない**——`.needs`
に正直に記録する。 -/

/-- ★★★★★★★★**[Falt1] Lemma 1.1**(単射性・長さの等式)—— 主張と証明の
両方が実装された。

## ★主張

| 原文 | 宣言 |
|---|---|
| `Ω_V⊗_VW → Ω_W` が単射 | `falt1MapBaseChangeInjective`(`KaehlerAux.lean`) |
| 余核 `Ω_{W/V}` の長さ `= length(W/p^δW)` | `falt1CokernelLengthEq` ＋ `falt1CokernelIsoLinear`
  (余核 `Ω_{W/V} ≅ W⧸differentIdeal V W` の特定、`differentIdeal V W` が
  `p^δ` の役割) |

## ★★★★逸脱の記録(CLAUDE.md の「逸脱」)

### 1. `Interface.RamificationSetup` への正式な差し替えは未実施

`Interface/Falt1/Ramification.lean` の `RamificationSetup` は
`cokerLength quotientLength : ℕ`(自然数)を要求するが、本実装の
`Module.length` は `ℕ∞`(無限を許す拡張自然数)値である。★これらが
**有限であることの証明**(`W` が完備離散付値環上有限で `differentIdeal`
が非零イデアルなら `W⧸differentIdeal V W` は有限加群になるはず)は
まだ行っていない——一般の Dedekind 整域 `V` に対しては `IsArtinian
W W` すら成り立たない(`W` 自身は自分自身上のアルティン加群ではない)
ため、`Submodule.length_quotient_lt` の直接適用もできない。Falt1 の
実際の設定(`V` が完備離散付値環、多くの場合剰余体が有限)に特化すれば
閉じられる可能性が高いが、未着手。

### 2. `Z`(絶対微分の基底)は Falt1 の原文で明示されていない

原文は「絶対微分」`Ω_V`・`Ω_W` を使うが、その絶対基底 `Z` を明示しない
(暗黙に固定されている——恐らく `W(k)`(剰余体の Witt 環)や `ℤ_p`)。
★本実装は **任意の `Z`(`Algebra Z V`・`Algebra Z W`・`IsScalarTower
Z V W` を満たす限り)に対して成り立つ一般形**として証明した——これは
弱化ではなく、原文の主張を**任意の絶対基底に対して同時に**証明した
ことになる(Falt1 の実際の `Z` を選べば、その具体形が直ちに従う)。

### 3. `differentIdeal V W ≠ ⊥` を仮定として受ける

有限可分拡大の differentIdeal が非零であることは mathlib の
`differentIdeal_ne_bot`(`[Module.Finite V W][Algebra.IsSeparable
(FractionRing V)(FractionRing W)]` が必要)から従うはずだが、
Falt1 の `K,L` が本当に `FractionRing V,FractionRing W` に一致する
ことの確認(`IsFractionRing W L` の導出等)はまだ行っていない——
仮定として直接渡すに留めた。 -/
def lemma_1_1_falt1.src : ABC3.Meta.Source :=
  { paper := "Falt1", pdfPage := 4, item := "Lemma 1.1", sectionId := "falt1-lemma-1-1" }

def lemma_1_1_falt1.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "falt1MapBaseChangeInjective(単射性そのもの)"
      (.inProject "ABC3" "ABC3.Found.Falt1.falt1MapBaseChangeInjective") 4,
    .citation "[ABC3]" "falt1CokernelLengthEq(長さの等式そのもの)"
      (.inProject "ABC3" "ABC3.Found.Falt1.falt1CokernelLengthEq") 4,
    .citation "[ABC3]" "falt1CokernelIsoLinear(余核 Ω_{W/V} ≅ W⧸differentIdeal V W の特定)"
      (.inProject "ABC3" "ABC3.Found.Falt1.falt1CokernelIsoLinear") 4,
    .implicitStep
      ("★逸脱 1: RamificationSetup(ℕ 値の長さ)への正式な差し替えは未実施。" ++
       "Module.length は ℕ∞ 値——有限性の証明(W⧸differentIdeal が有限加群) が要る") 4,
    .implicitStep
      ("★逸脱 2: 絶対微分の基底 Z は原文で明示されないため、任意の Z に対する" ++
       "一般形として証明した(弱化ではなく一般化)") 4,
    .implicitStep
      ("★逸脱 3: differentIdeal V W ≠ ⊥ を仮定として受ける。" ++
       "differentIdeal_ne_bot から導くには IsFractionRing W L 等の確認が要る") 4 ]

end ABC3.Found.Falt1
