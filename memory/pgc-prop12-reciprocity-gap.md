---
name: pgc-prop12-reciprocity-gap
description: 節目(5)(Lubin-Tate相互律の射影極限、sorry無しで完成)からProposition 1.2(Γ_K^ab≅(K^×)^∧)へ橋渡しするために追加で要るもの——未着手・未構築であることを確認(2026-09-05)
metadata:
  type: project
---

節目(5)(`memory/pgc-lubin-tate-existence-progress.md`——古典的
Lubin-Tate理論の射影極限による局所類体論の実現、`(i)principalUnits
QuotientEquivの自然性・(ii)reciprocityMapLimitの構成・(iii)その
全射性と核の特徴づけ`)が構成的にはすべて完成した。これを`/goal`の
残り8項目、特に**Proposition 1.2**(`lean/ABC3/Skeleton/PGC/Section1.
lean::residueCard_and_degree_recoverable`)へ接続できるかを検討した
結果を記録する。

## Proposition 1.2が実際に主張すること

`residueCard_and_degree_recoverable RD : (residueCardAndDegreeObject
RD).RecoverableFromAbsGal`——展開すると:

```
∀K K'(α:ContinuousMulEquiv K.absGal K'.absGal), (RD.card K,[K:ℚ_p]) =
  (RD.card K',[K':ℚ_p])
```

（`residueCardAndDegreeObject.transport := fun _ x => x`が恒等なため）
——**任意の位相群としての同型`Γ_K≅Γ_K'`に対し、剰余体の元の個数`q`と
絶対次数`[K:ℚ_p]`が一致する**、という**純粋に群論的な**主張。原文の
論拠(古典的な局所類体論`Γ_K^ab≅(K^×)^∧`とそこからの計数)は**証明の
一つの経路**であって、命題自体はこの経路を経由しなくても成り立てば
よい——ただし現時点で他の経路は見えていない。

## 節目(5)がカバーする部分と、まだ無い部分

節目(5)(`reciprocityMapLimit`)が実現するのは、**完全に分岐した
アーベル拡大`K_π=K(Λ_∞)`のGalois群**`Gal(K_π/K)`(=`Gal(K̄/K)`を
`ker(reciprocityMapLimit)`で割ったもの、`K_π`自体は中間体として
構成していない)が`𝒪_K^×`(の`CompatibleUnits`実現)と**同型**である
こと。古典的局所類体論の主張`Γ_K^ab≅(K^×)^∧`は、これに加えて:

1. **不分岐拡大の理論**: `K^ur`(`K`の最大不分岐拡大)・`Gal(K^ur/K)
   ≅Ẑ`(Frobenius元が生成する)。★2026-09-05実測:
   `.cache/mathlib-index.txt`・`.cache/decl-index.txt`(PGCフォルダ)
   のいずれにも、不分岐拡大・Frobenius元・`Gal(K^ur/K)`に相当する
   局所体理論の道具は**見当たらない**——`Mathlib/RingTheory/
   Frobenius.lean`にあるのは可換環論的な「算術Frobenius」
   (`IsArithFrobAt`)であって、局所体の不分岐拡大のGalois群という
   文脈のものではない。**この部分は完全に未着手・mathlib不在**。
2. **`K^×`の付値による分解**: `K^×≅𝒪_K^××π^ℤ`(位相群として、`π`を
   1つ固定すれば`v:K^×→ℤ`が分裂する)——これ自体は易しいはずだが
   本セッションでは未確認・未構築。
3. **上記(1)(2)と`Gal(K_π/K)≅𝒪_K^×`(節目(5))を組み合わせて
   `Γ_K^ab≅(K^×)^∧`(完備化込み)を結論する段**——`Γ_K^ab`が
   `Gal(K_π/K)`と`Gal(K^ur/K)`の(直積または拡大としての)組み合わせ
   になることを示す必要があり、これも未着手。
4. **`Γ_K^ab`(あるいは`(K^×)^∧`)の抽象的な群構造だけから`q`と
   `[K:ℚ_p]`を読み取る段**(Proposition 1.2の実際の主張に必要な
   最後の一歩): `𝒪_K^×≅μ_{q-1}×(1+𝔪_K)`(捩れ・pro-p部分への分解)
   から`q`(捩れの位数+1)と`[K:ℚ_p]`(pro-p部分の`ℤ_p`-階数−1)を
   取り出す議論——これも未着手。

## 結論・次に戻るときの判断

節目(5)は`Γ_K^ab≅(K^×)^∧`という古典的局所類体論の主張の**一部**
(完全分岐アーベル拡大の部分)を厳密に実現したが、Proposition 1.2
自体を閉じるには上記4点すべてが要り、**特に(1)不分岐拡大の理論は
mathlibに一切無く、本プロジェクトでも未着手**——これは節目(5)全体
(本セッション1日分の作業量)に匹敵するかそれ以上の規模になる見込み。

CLAUDE.mdの姿勢(工数の山を「壁」と呼ばない・既知数学のperson-years
は壁でなく道)に従い、これを新たな**スケルトンの候補**として記録
する——次にこの方向で戻る際は、`Interface/PGC/LocalFieldData.lean`
に「不分岐拡大」「Frobenius」を新しい`Interface`の自由なデータ(または
`Skeleton`の新規ファイル)として明示的に追加するところから始めるのが
筋(依存グラフを更新してから、層番号の低い葉から着手するという通常の
進め方)。

一方、Proposition 1.2以外の`/goal`残り項目(Corollary 1.3はProp 1.2に
依存・Proposition 2.1/2.2はp進対数+Verlagerungの経路で節目(5)とは
**別系統**・Corollary 3.1/3.3はHodge-Tate/uniformizing加群でこれも
別系統・Theorem 4.2は最終組み立て)は、節目(5)を直接には必要としない
可能性がある——特にProposition 2.1(`Section2.lean::prop_2_1`)は
`.needs`によればp進対数の3性質(準同型性・単射性・全射性)がすべて
本プロジェクトで**既に確立済み**(`Found/PGC/PadicLogMul.lean`・
`PadicLogInjective.lean`・`PadicLogSurjective.lean`)・Verlagerungも
mathlibに存在(`MonoidHom.transfer`)——境界外入力がほぼ解消済みで、
残るのは`prop_2_1.needs`の3番目`.implicitStep`(これらを組み合わせて
`RecoverableAsAddModule`を結論する段)のみという状態。**次にLubin-Tate
方面から離れて`/goal`側へ直接戻るなら、節目(5)より先にProposition 2.1
の方が近い可能性がある**——ただし原典の証明文が「独立したproofブロック
を持たず、直前の地の文がそのまま論拠になっている」ため、まず原典の
該当箇所(pGC p.4冒頭)を精読して論拠を正確に再構成する必要がある。
