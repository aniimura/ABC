# 1_Structured: 層B(②③一体)の対象——Proposition6.7・Definition6.4(i)・Definition6.1(iii)(iv)・Example4.4(i)(ii)

作成日: 2026-08-14(相談役指摘により事後的に整理——本セッションの読解調査[`3_Purified/.../proposition6.7-merge-welldefinedness.md`]は既に逐語抽出・行番号照合を伴っていたが、正式な1_Structured/2_LocatorMap文書として独立に整理されていなかった。工程を飛ばす既知のパターン[M7の3_Purified/4_CalibrationPlan飛ばし]の再発を避けるため、本文書で1_Structured相当を確定する)。

## 段落単位への分割(層B[②③一体]に関わる4段落)

### 段落A: Proposition6.7(全文、l.11765-11800)

「Then by replacing †DT by †DT⋇ [cf. Definition6.4,(i)], identifying the D-prime-strip †D≻ with the D-prime-strip †D0 via †φΘ±0 [cf. the discussion of Definition6.4,(i)] to form a D-prime-strip †D>, replacing the various +-full poly-morphisms that occur in †φΘ±± at the v∈Vgood by the corresponding full poly-morphisms, and replacing the various +-full poly-morphisms that occur in †φΘ±± at the v∈Vbad by the poly-morphisms described [via group-theoretic algorithms!] in Example4.4,(i),(ii), we obtain a functorial algorithm for constructing a [well-defined, up to a unique isomorphism!] D-Θ-bridge...」

**構成要素分解**(4要素、並列列挙構文「by [A], [B], [C], and [D], we obtain...」):
- A = T⋆置換(Definition6.4,(i)自体のT⋆構成、層A、他から独立、既に検証済み)
- B = ②(†D≻/†D0合体、†φΘ±0経由)
- C = ③-good(v∈VgoodでΘ±full→full置換)
- D = ③-bad(v∈VbadでΘ±full→Example4.4,(i),(ii)由来のpoly-morphismへ置換)

### 段落B: Definition6.4,(i)(全文、l.11296-11341、特にl.11326-11338)

「so each constituent D-prime-strip of †D|T| is only well-defined up to a positive automorphism, but this indeterminacy will not affect applications of this construction — cf. Propositions6.7; 6.8,(ii); 6.9,(i), below」

Proposition6.7への**直接の前方参照**を含む——B(②合体)がwell-definedである条件についての、Definition6.4(i)自身の予告的説明とみなせる。

### 段落C: Definition6.1,(iii)(iv)(l.10944-11002)

(iii): 「Aut+(†Dv)⊆Aut(†Dv)」(単一場所版、l.10944-10948)・「Aut+(†D)⊆Aut(†D)」(D-prime-strip全体版、l.10954-10957)——「positive automorphisms」の**唯一の定義箇所**。
(iv): 「+-full poly-isomorphism...obtained as the Aut+(†Dv)-...orbit of an isomorphism」(l.10964-10967)——置換前の†φΘ±0が定義上必ず全単射である根拠。

### 段落D: Example4.4,(i)(ii)(l.7113-7247、特にl.7183-7190)

「φΘ_vj...obtained by composing with...the natural surjection Πv↠Gv...from the evaluation sections labeled j」——③-bad(D要素)の実際の構成、一般に非単射。

## mathlib境界スキャン(軽量、2_LocatorMap相当)

- 「index-two subgroup」「kernel of a surjection onto {±1}」型の構造: mathlib`Subgroup.index_eq_two_iff`(既にセッション1〜2で発見・使用済み、`FINDING_POLY_MORPHISM_COLLAPSE_PROP67`解消に使用)。
- 「全単射性を性質として切り出す」設計(`Function.Bijective`): mathlib標準、既にtoy(`OpaqueProp67MergeStepToyCheck.lean`)で使用済み。
- 「Aut+(†D)のD-prime-strip全体版」(複数場所vにわたる自己同型の直積構造): mathlib`MulAut`・`Pi`型の組み合わせで表現可能と推測(既存`AutPlusDv`の単一場所版パターンを`∀ v`へ拡張する形)——次段(層B型設計本体)で具体的に確認が必要。
- 新規axiom化が必要な境界外候補: なし(現時点で発見された全ての構成要素は既存mathlib型+既存プロジェクト資産で表現可能と見込まれる)。

## 2_LocatorMap: 相互参照の解決

| 用語 | 定義箇所 | 再利用箇所 |
|---|---|---|
| positive automorphism / Aut+(†D) | Definition6.1,(iii)(l.10944-10957) | Definition6.4,(i)(l.11336-11338、引用符なし再使用、citation tagなし) |
| +-full poly-isomorphism | Definition6.1,(iv)(l.10964-10967) | Proposition6.7(l.11773-11780、置換前の†φΘ±0) |
| evaluation sections / Example4.4由来のpoly-morphism | Example4.4,(i),(ii)(l.7113-7247) | Proposition6.7(l.11779-11783、Vbad置換後の†φΘ±0) |
| †φΘ±0 | Definition6.4,(i)の議論部分(「cf. the discussion of Definition6.4,(i)」、Proposition6.7自身が参照) | Proposition6.7の②(合体)本体 |

## 引き渡し(2_LocatorMap→3_Purified、既に実施済み)

上記の全構成要素は`3_Purified/Inter-universal Teichmuller Theory I/proposition6.7-merge-welldefinedness.md`で既に精製済み——本文書はその1_Structured/2_LocatorMap相当を事後的に確定するものであり、内容の重複でなく工程記録の補完。
