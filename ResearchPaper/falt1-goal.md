# [Falt1] トラックのゴール(2026-09-04 設定・確定、Chapter I は 2026-09-04 完了)

対象: G. Faltings, *p-Adic Hodge Theory*, JAMS Vol.1 No.1 (1988) pp.255-299
(物理46ページ、`papers.json`に短縮タグ`Falt1`で登記済み、`0_Source`にPDF/txtあり)。

## 0. ゴール(現在地)

> **Falt1 Chapter I §1 2/2, §2 4/4, §3 2/2, §4 5/5 —— Chapter I 13/13 達成(2026-09-04)**
>
> 残り: Chapter II §1 0/2, §2 0/4, §3 0/2,
> Chapter III §1 0/3, §2 0/2, §3 0/3, §4 0/1 —— 合計 0/17。
> 全体では 13/30。

★★**正直な注記(過大申告を避けるため)**: Chapter I 13項目はすべて
`Interface`(posit、証明していない)+`Skeleton`(.src・.nonvacuous 付き)
+ 260dpi目視確認+ HTML逐語アンカーまで揃い、`node tools/check.mjs`・
`node tools/graph.mjs --owner Falt1`(条なし13・項目13・sorry0)を通る。
ただし:
- §4.4(Theorem 4.4)は[LocProP]形式化時に作った**(i)(iv)のみの縮小形**
  のまま(原文は(i)-(vii)まである)。`.needs`に理由を明記済み。
- §4.2・§4.5は原文が持つ添字(i・M・R1,R2 の分解等)に関する**全称量化を
  代表的な1インスタンスに簡略化**して posit した(`Interface/Falt1/
  GaloisCohomologyCompute.lean` 冒頭に記録)。
- 「almost mathematics」の核心である almost étale covering・m-torsion
  等は mathlib に対応物が無いため全項目が posit ベース(`sorry` は無いが
  `Found`(証明)ではない)。数学的に**証明した**わけではない。

★★★**2026-09-04、`/goal` で「Chapter I を Found(証明)にする」という
より厳しい条件が設定された——その進捗と実現可能性の評価を追記する。**

## 0.1 `/goal Falt1 Chapter I Found` の進捗(2026-09-04)

**Lemma 1.1**(§1 1/2): ★★★★★**完成**。`Found/Falt1/KaehlerAux.lean`・
`Found/Falt1/Lemma11.lean` として数学的内容(単射性・余核の長さの等式)
を**完全に sorry無しで証明した**(`falt1MapBaseChangeInjective`・
`falt1CokernelLengthEq`・`lemma_1_1_falt1`)。途中で作った再利用可能な
一般補題(mathlib に無かったもの): `subsingleton_H1Cotangent_self`
(`FormallySmooth R S ⟹ Subsingleton(H1Cotangent R S)`)・
`polynomialKaehlerSplit`(`Ω_{V[T]/Z} ≅ (V[T]⊗_VΩ_V)×V[T]` の直和分解)・
`mapBaseChange_injective_transport`(`mapBaseChange` の単射性を代数同型
に沿って輸送)。`differentIdeal_ne_bot` の全条件を `falt1_
differentIdeal_ne_bot` で導出——`differentIdeal V W ≠ ⊥` は仮定では
なく定理になった。★★★(2026-09-04 完成)`Interface.RamificationSetup`
への正式な差し替えも完了した(`falt1RamificationSetup`)——`Module.
length` の `ℕ∞→ℕ` 有限性は当初の見積り(CRT+局所化不変性)より単純に
出た(`Ring.DimensionLEOne`+`Ideal.krullDimLE_zero_quotient_iff_
forall_minimalPrimes_isMaximal`経由)。`RamificationSetup` の全
フィールドを本物のデータで埋め、`lemma_1_1 (falt1RamificationSetup
w hint hadjoin hw)` が型検査を通ることを確認した——Skeleton の posit
(`RamificationSetup.example`)を初めて実データに置き換えた例。
残る記録事項は「絶対微分の基底 `Z`」を一般化(`Z:=ℤ` を正準に選択)
として証明したことのみ(弱化ではない、`.needs` に記録)。

**Theorem 1.2**(§1 2/2): ★★**未着手。260dpi で物理p.5(印字p.258)を
直接読み(pdftotext の OCR は数式記号を激しく壊すため、上記 .txt 行
223-235 は信頼できない——実際に画像で読み直した、2026-09-04)、
実現可能性を評価した結果、Lemma 1.1 より明確に難しいことが分かった。**

- 設定: `V=V_0⊂V_1⊂V_2⊂...` という拡大の列で、`Ω_{V_n/V_{n-1}}`
  (★添字は `V_{n-1}` であって `V` ではない——OCR 行 223-235 の誤読を
  訂正)が `(V_n/pV_n)^{d+1}` を商に持つという条件を満たすもの(典型例:
  `V = W(k)(T_1,...,T_d)` の完備化、`V_n` は `T_i` の `p^n` 乗根と
  `p^{n+1}` 乗根の1のべき根を添加した拡大——perfectoid 系列の原型)。
- 主張: `W_n :=` `V_n⊗_VW` の正規化として、`W_n/V_n` の different
  `δ_n` が `n→∞` で `0` に収束する。
- 証明(原文どおり、正確に読み直した): 塔 `V_n→V_{n+1}→W_{n+1}`
  (と `W_n→W_{n+1}`)から `Ω_{W_n/V_n}⊗_{W_n}W_{n+1} →
  Ω_{W_{n+1}/V_n} → Ω_{W_{n+1}/V_{n+1}} → 0`(完全、これ自体は既存の
  `KaehlerDifferential.exact_mapBaseChange_map` で足りる)。
  **★★新たに要る部品**: 「`Ω_{W_{n+1}/V_n}` が `W_{n+1}/p^αW_{n+1}`
  型の加群 `d+1` 個の直和である」という構造的事実——これは Lemma 1.1 を
  `V_n→V_{n+1}` の**各生成元ごとに**適用して得られるはずだが、原文は
  この構成を「明らか」として省略しており、**単一生成元版の Lemma 1.1
  (今回証明した形)をそのまま複数生成元へ一般化する新規の仕事**になる。
  その後、`Module.length` の完全列に対する加法性・不等式(カーネル・
  コカーネルの長さ評価)を連鎖させ、最終的に `δ_n-δ_{n+1} ≥
  β-(d+1)(δ_n-δ_{n+1})`(`β=min{1,δ_n/(d+1)}`)という漸化不等式から
  `δ_n→0` を導く。

★★★**この構造から、Theorem 1.2 は Lemma 1.1 と同等かそれ以上の
新規インフラ(複数生成元版 Lemma 1.1・`Module.length` の不等式操作の
連鎖・「非常に分岐した」列の定義そのものの形式化)を要求する**——
Lemma 1.1 1件の形式化がこのセッションの大半(数十回のラウンド)を
要したことを踏まえると、Theorem 1.2 は**単独でさらに同程度以上の
作業**になる見込み。★急いで低品質な形式化を試みるより、正確に
読み直した上でこう記録する方が誠実だと判断した。

★★★2026-09-04、**「複数生成元版 Lemma 1.1」の核心(2生成元の場合)が
完成した**(`Found/Falt1/KaehlerAux.lean` の `pushoutKaehlerSplit`)。
`Algebra.IsPushout R C D B`(`B` が `C`・`D` の `R` 上の pushout=
テンソル積)のとき、`Ω_{B/R} ≃ₗ[B] (B⊗_CΩ_{C/R}) × (B⊗_DΩ_{D/R})`
(`mapBaseChange R C B` が単射という条件つき)——`polynomialKaehlerSplit`
と全く同じ「`Function.Exact.splitSurjectiveEquiv` + 明示的セクション」
のパターンで構成できた。**鍵になった発見**: `tensorKaehlerEquivBase`
(片方の脚基点)ではなく **`tensorKaehlerEquiv`**(もう1つの、`B⊗_D
Ω_{D/R} ≃ₗ[B] Ω_{B/C}` という形の mathlib 既製品)の方が正しい道具
だった——`map R C B B ∘ mapBaseChange R D B = tensorKaehlerEquiv
R C D B` という自然性(`tensorKaehlerEquiv_eq_map_mapBaseChange`、
新規に証明)が、`polynomialEquiv` の代わりのセクションを直接与えた。

★残る作業(次のセッションはここから):
1. **`d+1` 個への一般化**(今は2個のみ)——`pushoutKaehlerSplit` を
   `d+1` 回の pushout の反復に適用する帰納法(構造自体は難しくないはず)。
2. ★★★★**2026-09-04、`mapBaseChange R C B` の単射性の条件——
   完成した**(`mapBaseChange_injective_adjoinRoot_tensor`)。下の
   「3つ目(AdjoinRoot基底変換)」の経路が閉じた: `tensorPolynomialAlgEquiv
   (C⊗[R]Polynomial R ≃ₐ[C] Polynomial C)` → `adjoinRootTensorEquiv
   (C⊗[R]AdjoinRoot g ≃ₐ[C] AdjoinRoot(g.map(algebraMap R C)))` →
   `mapBaseChange_injective_of_nzd` + `mapBaseChange_injective_transport`
   で結ぶ。**`B = C⊗[R]AdjoinRoot g` の場合、`pushoutKaehlerSplit` の
   `hinj` は Lemma 1.1 とまったく同じ形の非零因子条件 1 つ
   (`g.map(algebraMap R C)` の微分が `AdjoinRoot(g.map(algebraMap R C))`
   で非零因子)に帰着する。** lake build 済み・sorry 皆無確認済み
   (commit `f45fa30b`)。残る技術的な下ごしらえ: `PUnit` の宇宙を
   `Type`(宇宙0)に固定する必要があった(`R C` を `Type*` にした瞬間
   宇宙メタ変数 `?u` が解決不能になる実測——`tools/lean-idioms.md` へ
   追記の価値あり)。
   ★★2026-09-04(この項目は上で解消済み、経緯として残す)経路を検討・
   部分的に実測した: `Algebra.
   Extension.cotangentComplex_injective_iff P [FormallySmooth R P.Ring]
   : Function.Injective P.cotangentComplex ↔ Subsingleton(H1Cotangent
   R A)`(mathlib)を `P := Extension.ofSurjective(algebraMap
   (Polynomial V)(AdjoinRoot g))` に適用すれば `[FormallySmooth V
   (Polynomial V)]`(既にある)から `Subsingleton(H1Cotangent V
   (AdjoinRoot g)) ↔ Function.Injective P.cotangentComplex` が出る
   はず——問題は **`P.cotangentComplex` が既存の `kerCotangentToTensor`
   と型レベルで(`rfl` でも)一致しないこと**(`P.Cotangent` は定義上
   `P.ker.Cotangent` で `P.ker` も定義上 `RingHom.ker(algebraMap
   P.Ring S)` のはずだが、`Extension.ofSurjective` の内部実装が
   `algebraMap` を直接使わないため、期待通りに簡約されない)。両者を
   **関数として**(`⇑` レベルで)結びつける橋渡しが未解決——`Extension`
   の内部実装(`RingTheory/Extension/Cotangent/Basic.lean` 等)を
   もう少し詳しく調べる必要がある。

   ★★★2026-09-04、**別の(より具体的な)経路も検討した**: `B = C⊗_R
   D`(`D=AdjoinRoot g`)は `C⊗_RPoly R ≅ Poly C`(`MvPolynomial.
   algebraTensorAlgEquiv`、単変数の場合に相当)+ テンソル積が商と
   可換という事実を経由すれば **`B ≅ AdjoinRoot(g.map(algebraMap R C))`
   ——つまり `B` は `C` 上の monogenic 拡大そのもの**になるはず。
   これが成り立てば `mapBaseChange R C B` の単射性は
   `mapBaseChange_injective_of_nzd`(`V:=C`・`f:=g.map(algebraMap R C)`)
   を**そのまま再利用できる**(H1Cotangent 経由の抽象論より遥かに
   具体的)——ただし `C⊗_RPoly R ≅ Poly C` の同型と「テンソル積は
   商と可換」を貼り合わせて `B ≅ AdjoinRoot(...)` を確立する作業
   自体は未着手(`MvPolynomial.algebraTensorAlgEquiv` を単変数
   `Polynomial` へ specialize する下ごしらえが要る)。**3つの経路
   (H1Cotangent+平坦底変換・Extension.cotangentComplex 橋渡し・
   AdjoinRoot 基底変換の同一視)を検討し、いずれも即座には閉じ
   なかった**——この時点ではここで打ち切った。
   →★★★★続きは同日内、上の「完成した」を見よ:3つ目
   (AdjoinRoot 基底変換)を実際に最後まで実行し、閉じた。
★★★★★2026-09-04、**Theorem 1.2 の証明全文(物理p.5=印字p.258)を
260dpi で precisely 読み直した**(OCR ではなく目視)。以下、逐語に近い
形で構造を記録する(今後の形式化がここから始められるように)。

> Proof. Replacing V by some V_n we may assume that all W_n are integral
> domains. Look at the sequence of maps Ω_{W_n/V_n}⊗_{W_n}W_{n+1} →
> Ω_{W_{n+1}/V_n} → Ω_{W_{n+1}/V_{n+1}}. The kernel of the second map
> contains Ω_{V_{n+1}/V_n}⊗_{V_{n+1}}W_{n+1}, which has (W_{n+1}/pW_{n+1})^{d+1}
> as quotient. As Ω_{V_{n+1}/V_n}⊗_{V_{n+1}}W_{n+1} is the direct sum of
> d+1 modules of the form W_{n+1}/p^αW_{n+1}, the kernel of the second
> map contains the kernel of multiplication by p on Ω_{W_n/V_n}⊗_{W_n}W_{n+1},
> and hence the composition of the two maps annihilates the kernel by
> p-multiplication on Ω_{W_n/V_n}⊗_{W_n}W_{n+1}. If we denote the different
> of W_n over V_n by δ_n, this kernel has length at least equal to that
> of W_{n+1}/p^βW_{n+1}, where β = min{1, δ_n/(d+1)}. Also it is clear
> that p^{δ_n-δ_{n+1}} annihilates W_{n+1}/(W_n⊗_{V_n}V_{n+1}), and hence
> the cokernel of the composition of the two maps is annihilated by
> p^{δ_n-δ_{n+1}}. So its length is at most equal to that of W_{n+1}
> divided by the (d+1)st power of this. We derive that
> δ_n-δ_{n+1} ≥ β-(d+1)(δ_n-δ_{n+1}). So if δ_n ≥ d+1, then
> δ_{n+1} ≤ δ_n-1/(d+2), and otherwise δ_{n+1} ≤ (1-1/((d+1)(d+2)))δ_n.
> In any case δ_n → 0 for n→∞. □

直後の「典型例」段落: `V = W(k)(T_1,...,T_d)` の完備化(`W(k)` は完全体
`k` の Witt ベクトル)、`V_n` は「`T_i` の `p^n` 乗根 + `1` の `p^{n+1}`
乗根」で生成される拡大——つまり **`V_{n+1}/V_n` は `d+1` 個の
単生成拡大(各 `T_i` の追加の `p` 乗根 1 個 + 冪根 1 個)の
**同時**添加**(逐次ではない)。これは `pushoutKaehlerSplit` の
`d+1` 項版が直接与える「`Ω_{V_{n+1}/V_n}` は `d+1` 個の巡回加群の
直和」という一文と**厳密に一致する**——★上の「3つの経路」で
確立した `mapBaseChange_injective_adjoinRoot_tensor` と合わせて、
`pushoutKaehlerSplit` は当初の想定通り Theorem 1.2 の証明の**まさに
その箇所**(第3文)で使われるべき道具であることが、憶測ではなく
原文の記述と照合して確認できた。

★以上を踏まえ、残りの形式化作業を(既存の(3)(4)を置き換えて)
以下のように精緻化する:

3a. **`pushoutKaehlerSplit` の `d+1` 項化**——`Fin (d+1)` 添字族
    `D : Fin (d+1) → Type*`(各 `R`-代数、`Algebra.IsPushout` の
    反復合成で `B := D 0 ⊗_R D 1 ⊗_R ⋯ ⊗_R D d` を構成)に対し
    `Ω_{B/R} ≃ₗ[B] ⊕ᵢ (B ⊗_{D i} Ω_{D i/R})`。★この段階では純粋に
    Kähler微分の代数——`V_n`族そのものへの依存は無い、汎用インフラ
    として先に作れる。

    ★★★★2026-09-04、**`d=2`(3項)の場合を最後まで実行し、帰納の
    1ステップが機械的に閉じることを実証した**(`pushoutKaehlerSplit3`、
    commit `d8341aab`)。`pushoutKaehlerSplit` を2回・
    `TensorProduct.AlgebraTensorModule.cancelBaseChange`(「反復底変換
    =直接底変換」)・`TensorProduct.prodRight`(テンソル積と直積の
    可換性)の組み合わせで構成できた——`d+1` 項への一般化もこの2つの
    道具の反復のみで閉じる見込みが実証された。`Algebra C0 B`・
    `Algebra C1 B` とその `IsScalarTower` は(`pushoutKaehlerSplit`
    自身と同じ流儀で)**仮定として明示的に渡す**必要がある——`B1`・`B`
    を具体的な `TensorProduct` として書いた特殊化版は、宣言の**型
    そのもの**が `Algebra.TensorProduct.rightAlgebra`(`abbrev`、
    tools/lean-idioms.md #23)を要求するため別途書けないことも実測した
    (使用箇所ごとに `letI := Algebra.TensorProduct.rightAlgebra` を
    先に置いてから直接呼ぶ必要がある)。
    次のセッションへ: `Fin (d+1)` 一般の帰納は、**型族を再帰的に
    構成する**部分(`B : ℕ → Type*` を `TensorProduct` で再帰定義し、
    各段の `CommRing`/`Algebra` インスタンスを同時に構成する)が
    Lean4 では地雷("recursive Type-valued def with dependent
    instances")——`structure RAlg (R) where carrier; [ring];[alg]`
    のような Σ 束ねを検討すること(未着手)。

    ★★2026-09-04、この一般化への着手を**意図的に保留した**——下の 3b
    の難度を精密に再評価した結果、`pushoutKaehlerSplit` の `d+1` 項化
    (純粋な Lean エンジニアリング)を仕上げても、3b(独立の length
    計算、Lemma 1.1 同等の数学的深さ)が閉じない限り Theorem 1.2 全体
    には接続しない、と判断したため。3a を一般化する労力の投資対効果は
    3b の見通しが立ってから再評価する。
3b. **`Ω_{V_{n+1}/V_n}⊗W_{n+1}` が `d+1` 個の巡回加群の直和という
    事実そのものを `V_n` 族の定義に接続**(3a の適用)。

    ★★★★★2026-09-04、**この項目の難度を精密に再評価した**(原文
    証明の中段落を再度 260dpi で読み直し、独立に再構成を試みた結果)。
    `pushoutKaehlerSplit3` が与えるのは「直和分解の**存在**」だが、
    `delta_tendsto_zero` が要求する `hrec`(漸化不等式)を実際に導く
    には、原文の中段落——「第1の写像(`Ω_{Wₙ/Vₙ}⊗Wₙ₊₁→Ωₙ₊₁/Vₙ`、
    タワー `Vₙ→Wₙ→Wₙ₊₁` の base change)と第2の写像(`Ωₙ₊₁/Vₙ→
    Ωₙ₊₁/Vₙ₊₁`、`KaehlerDifferential.exact_mapBaseChange_map` により
    その kernel はタワー `Vₙ→Vₙ₊₁→Wₙ₊₁` の base change の**像**と一致)
    の**合成**の kernel・cokernel を、`δₙ`(`Wₙ/Vₙ` の different の
    長さ)を使って評価する」という**独立の length 計算**が要る。
    特に「`p^{δₙ-δₙ₊₁}` が `Wₙ₊₁/(Wₙ⊗_{Vₙ}Vₙ₊₁)` を零化する」という
    一文は、`Wₙ⊗_{Vₙ}Vₙ₊₁` と(その正規化である)`Wₙ₊₁` の差を
    differentで測る、半局所環の conductor 理論に相当する事実で、
    **Lemma 1.1 自体と同等かそれ以上の深さの独立した数学的主張**
    であり、`pushoutKaehlerSplit3` を「適用するだけ」では閉じない
    ——単なる工学的な貼り合わせ作業ではなく、**新しい数学的補題
    (半局所 Dedekind 的環の塔での conductor/different の劣加法性)
    を要する**、という結論に達した。★安易な誤形式化を避けるため、
    この評価が固まるまで 3b への Lean 実装は保留する。

    ★★★★★★★2026-09-04、**上記の等式を実際に Lean で構築・証明した**
    (`differentIdeal_tower_diamond`、`Found/Falt1/KaehlerAux.lean`、
    commit `937de066`、sorry無し)——`Algebra.IsSeparable` の引数を
    `FractionRing.liftAlgebra` に明示的に固定しないと2回目の適用で
    instance 探索が食い違うことが判明し(tools/lean-idioms.md #23と
    同種)、それを回避して完成させた。下記の分析(未知数
    `differentIdeal Wₙ Wₙ₊₁` の評価が核心の困難として残る)は変わらず
    有効——これは Theorem 1.2 3b の**代数的骨格**として独立した
    検証済みの一歩である。

    ★★★★★★2026-09-04(発見時のメモ、上で実装完了): mathlib に
    **`differentIdeal` の塔の公式**が既にある
    (`differentIdeal_eq_differentIdeal_mul_differentIdeal`、
    `RingTheory/DedekindDomain/Different.lean`)——`A→B→C` の塔に対し
    `differentIdeal A C = differentIdeal B C * (differentIdeal A B).
    map(algebraMap B C)`(`A` 整閉整域・`B,C` Dedekind・有限・捩れ無し・
    `Algebra.IsSeparable (FractionRing A)(FractionRing C)` が仮定)。
    これを `Vₙ→Vₙ₊₁→Wₙ₊₁` と `Vₙ→Wₙ→Wₙ₊₁` の2通りの塔に適用すると
    `δₙ₊₁ · (differentIdeal Vₙ Vₙ₊₁).map(...) = (differentIdeal Wₙ
    Wₙ₊₁) · δₙ.map(...)` という**厳密な等式**(原文の不等式より強い形)
    が直接得られる——ただし残る未知数 `differentIdeal Wₙ Wₙ₊₁`
    (`Wₙ⊗_{Vₙ}Vₙ₊₁` が正規化前の状態からどれだけ `Wₙ₊₁` と異なるかを
    測る量)を評価する必要は結局残り、これが原文の「kernel・cokernel
    の length 評価」の**まさにその部分**に相当する——つまりこの公式は
    3b の**別の切り口**を与えるが、核心の困難(非正規性の評価)は
    避けられない。ただし「厳密な等式」から出発できる分、不等式の
    往復より見通しが良くなる可能性がある。次のセッションでの
    検討候補として記録する(まだ Lean での検証・使用はしていない)。
★★★★★★★★2026-09-04、**決定的な発見**: `ResearchPaper/0_Source/Brinon
Conrad - CMI Summer School Notes on p-adic Hodge Theory.txt` の
Exercise 13.7.4 が、Faltings のこの証明を**教科書レベルの精密さで
step-by-step に書き下している**(Brinon-Conrad はスタンフォード/
プリンストンで定番の CMI サマースクールノート)。原文(Faltings 自身
の1988年の論文)の圧縮された記述より遥かに読みやすく、以下の疑問点
が解消した:

- (1) `Ω¹_{B1/A}` が `d+1` 個の元で生成される、という事実は
  **「second fundamental exact sequence + Nakayama」で出す**(具体的な
  手法が明示されている——`pushoutKaehlerSplit` 系ではなく、もっと
  直接的な「単項生成の完全列」の議論)。
- (2) **`length_B(Ω¹_{B/A}) = length_B(B/p^{δ_{B/A}}B)`**(`δ_{B/A}`
  は different `D_{B/A} = p^{δ_{B/A}}B` の指数)——これは
  `Found/Falt1/Lemma11.lean` の `falt1CokernelLengthEq` と**本質的に
  同じ計算**(Lemma 1.1 そのもの)であることが確認できた。
- (4) `ker(b) ⊇ ker(p 倍写像)` の正確な意味と根拠が判明:
  「(1)(=d+1個生成)+ elementary divisor theorem」から出る——
  `Ω¹_{An+1/An}⊗Bn+1` が(有限生成 torsion 加群の構造定理により)
  `d+1` 個の巡回加群の直和に分解でき、各巡回加群の位数が `p^{何か}`
  なので、"合成 `b∘a` の kernel は `p` 倍写像の kernel を含む"
  という主張の正確な機構が判明した(★私が「よく分からない」と記録
  していた箇所)。
- (5) **鍵となる未解決の箇所の出典が判明**: 「using the definition of
  the discriminant, show `p^{δn-δn+1}·(Bn⊗_{An}An+1) ⊆ Bn+1`」——
  これは**判別式(discriminant)の塔での乗法性 + 非極大な order の
  conductor が判別式の比較から評価できる、という古典的整数論**
  (conductor-discriminant の関係)を使う、と明記されている。
  `differentIdeal_tower_diamond`(different の塔の公式)ではなく
  **discriminant の塔の公式**が本来の道具だった可能性が高い——
  次のセッションでは mathlib の `Algebra.discr`/`NumberField.discr`
  周りの資産を調査すること(未着手)。

★結論(更新): 核心の困難(非正規性の評価)は依然として**独立の
古典的整数論(conductor-discriminant の関係)を要する**という
これまでの評価は変わらないが、**その具体的な形(discriminant の
塔の乗法性)が判明した**ことで、次にどの mathlib 資産
(`Algebra.discr` 系)を調査すべきかが明確になった。これは
「思ったより深い」から「深いが具体的に何を調べればよいかが分かった」
への前進である。

★★★★★★★★★2026-09-04、続けて mathlib を調査し、**まさにこの
「conductor via discriminant/微分」の道具そのもの**を発見した:
`conductor_mul_differentIdeal`(`RingTheory/DedekindDomain/
Different.lean`)—— `(conductor A x) * differentIdeal A B =
Ideal.span {aeval x (derivative (minpoly A x))}`(`x : B` が
`Frac(B)/Frac(A)` を生成する場合)。これは Lemma 1.1 で使った
`differentIdeal_eq_span_derivative`(`conductor = (1)` の特殊形、
`B` が `A` 上 monogenic な場合)の**一般化**であり、`conductor A x`
が `A[x]`(= 私の `T_{n+1} = Wₙ⊗_{Vₙ}Vₙ₊₁` に相当、`x` が「典型例」
の生成元(`p`乗根等)なら monogenic)が `B`(= `Wₙ₊₁`)からどれだけ
離れているかを**まさに測る量**。つまり:

`differentIdeal Wₙ Wₙ₊₁ = Ideal.span{aeval x (deriv(minpoly Wₙ x))}
/ (conductor Wₙ x)`(概念的に——実際には `conductor * different =
span{...}` の形で使う)

という式で、目標の `Jₙ := differentIdeal Wₙ Wₙ₊₁` が
`conductor Wₙ x` 経由で書ける。★「典型例」(`V_{n+1}/V_n` が
`p`乗根の同時添加)なら `x` の具体的な選び方も明確。

★★★★★★★★★★2026-09-04、**実際に組み合わせて完成させた**
(`cancel_conductor_delta`、`Found/Falt1/KaehlerAux.lean`、commit
`bf0a3ea7`、sorry無し)。`differentIdeal_tower_diamond` と
`conductor_mul_differentIdeal` の2つの等式を掛け合わせるだけで `Jₙ`
が消え、**`conductor(Wₙ,x) · δₙ₊₁ = δₙ.map(...)` というきれいな
関係式**が(「`conductor Wₙ x` と `differentIdeal Vₙ Vₙ₊₁` の
base change が一致する」という `hspan_eq` 仮定のもとで)出ることを
確認した。Dedekind整域の「0でないイデアルは消去可能」という性質
だけで閉じる、純粋に代数的な議論。

★残る作業: (a) `hspan_eq`(3c=「典型例」で `x` の最小多項式が base
change と両立することの確認)、(b) `conductor(Wₙ,x)` 自体の下限
評価(`p^{何か}` を含むことの証明、Brinon-Conrad step (5) 後半に
対応)、(c) この関係式(掛け算の形)から `delta_tendsto_zero` が
要求する `hrec`(引き算の不等式の形)への変換。

★★★★★★★★★★★2026-09-04、(c) に実際に着手し、**「(a)(b)より軽い」
という見立てが誤りだったと判明した**。`length(R/(I·J)) =
length(R/I) + length(R/J)`(Dedekind整域、`I,J` 0でない、必ずしも
coprimeでない一般の場合)という長さの加法性が要るが、これは
mathlib に無い(coprimeの場合の `Ideal.quotientMulEquivQuotientProd`
はあるが一般の場合は無い)。証明の筋は分かる——完全列
`I/(IJ) → R/(IJ) → R/I → 0`(`Module.length_eq_add_of_exact` が
使える)+ `I/(IJ) ≅ R/J`(可逆イデアル `I` によるテンソルが
length を変えないという事実、局所化して「可逆イデアルは各素点で
自明」に帰着させる議論)——だが後半の同型は**局所化と length の
両立性**という、これも mathlib に無さそうな一般論を要する。
★このセッションで繰り返し遭遇したパターン(「見た目は軽いはずが
掘ると独立の一般補題に当たる」)がここでも再現した——安易に
「軽い」と決めつけず、正直に記録する。次のセッションへ:
`length(R/(I·J)) = length(R/I)+length(R/J)` を独立の補題として
先に確立すること(Dedekind整域論のごく基本的な事実のはずだが、
mathlib での正確な組み立て方は未確認)。

    ★★2026-09-04、この長さの加法性を埋める手がかりを追加調査した。
    `Ideal.relNorm : Ideal S →*₀ Ideal R`(`RingTheory/Ideal/Norm/
    RelNorm.lean`、**任意の Dedekind 底環 `R`** に対する相対ノルム、
    `→*₀` として bundled——`ℤ` 限定の `Ideal.absNorm` と違い一般の塔
    に使える)は型として乗法的性なので `relNorm(IJ)=relNorm(I)*
    relNorm(J)` はタダで手に入る。ただしこれは `S` 上の長さではなく
    `R` の**イデアルへ**写す不変量なので、`length_S(S/I)` を直接
    与えるには別途「`length_S(S/I) = length_R(R/relNorm I)`」的な
    比較公式が要り、それ自体が `R` 側で**同じ加法性**を要求する
    (循環)——この経路は今回は閉じなかった。次のセッションへ:
    「単一の素イデアル `P` に対し `length_S(S/P^n) = n·length_S(S/P)`」
    を最初の一段として直接攻める(`TensorProduct.tensorQuotEquivQuotSMul`
    で `P/P^{n+1} ≅ S/P` を示し帰納する筋、`P` が可逆という事実を
    経由)方が、`relNorm`/`absNorm` 経由より見通しが良いかもしれない。

    ★★★★★★★★★★★★2026-09-04、**核心の道具を発見した**:
    `RingTheory/PicardGroup.lean` の `Module.Invertible`(可逆加群、
    Picard群論のための typeclass、`contractLeft` の全単射性で定義)に、
    まさに必要な「テンソルは完全性を保つ」補題群がある:
    `Module.Invertible.rTensor_injective_iff`・`_surjective_iff`・
    `_bijective_iff`(`[Module.Invertible R M]` のとき
    `f.rTensor M` の単射性・全射性・全単射性は `f` 自身のそれと同値)。
    これは「不変(=可逆)加群でテンソルすると完全列が完全列のまま」
    という、`Module.length_eq_add_of_exact` と組み合わせれば
    `length(M⊗N)=length(N)` を**帰納法で**導出できるはずの、
    まさに探していた一般論——mathlib に無いと思っていたが、
    「Picard群」という一見無関係な名前のファイルに実在した
    (`measure-mathlib-before-skeleton` の教訓が今回も再現: 概念名
    ではなく数学的内容で探すべきだった)。

    ★残るブリッジ: `IsDedekindDomain S`(の nonzero イデアル `I`)から
    `Module.Invertible S ↥I` インスタンスを得る経路——
    `isDedekindDomain_iff_isDedekindDomainInv`(`RingTheory/
    DedekindDomain/Ideal/Basic.lean`)経由で `I * I⁻¹ = 1`
    (`FractionalIdeal` として)は出るが、これを `Module.Invertible`
    の `contractLeft` 全単射性の形へ変換する配線はまだ未着手。

    ★★★★★★★★★★★★★2026-09-04、**この一般化の代わりに、単項イデアルの
    場合に限定して完成させた**(`length_quotient_span_singleton_mul`、
    `Found/Falt1/KaehlerAux.lean`、commit `25221786`、sorry無し)。
    `Module.Invertible`(可逆加群・局所化)の一般論は一切不要で、
    既存の道具(`Module.length_eq_add_of_exact`・`Submodule.
    quotientQuotientEquivQuotient`・`LinearMap.quotKerEquivRange`)
    の組み合わせだけで、`a` が非零因子のとき
    `length(R/((a)·J)) = length(R/(a)) + length(R/J)` を証明できた。
    多くの実際の場面(`differentIdeal_eq_span_derivative`・
    `conductor_mul_differentIdeal` 等、単一の生成元の微分で書ける
    場合)ではイデアルが単項なので、3b(c) の核心 gap をこの特殊形で
    埋められる見込みが高い。★これで3b(c)は実質的に解決した——
    残るのは `cancel_conductor_delta` の各イデアル(`conductor Wₙ x`
    等)が実際に単項であることの確認(3a と同様、「典型例」の具体形
    次第)。

    ★★★★★★★★★★★★★★2026-09-04、残る (a)(b) の依存関係を精査した。
    (a)(`minpoly Wₙ x = (minpoly Vₙ の生成元).map(...)`という
    base change 整合性)は `minpoly.eq_of_irreducible_of_monic`
    (mathlib、体上でのみ成立)を経由するのが自然だが、これには
    「`f.map(algebraMap Vₙ Wₙ)` が `Frac(Wₙ)` 上既約のまま」という
    非自明な仮定が要り、これは**抽象的に証明できる一般命題ではなく、
    具体的な `Vₙ`・`Wₙ` の構成(3c)に依存する事実**である(一般には
    偽になりうる——`Wₙ` が `Vₙ` の拡大の情報を「知りすぎている」と
    既約性が壊れる)。(b)(conductor の下限評価)も同様に、具体的な
    ramification の値(`p` の分岐指数等)を要する。
    ★★結論: **(a)(b) は 3b の中の独立した補題としてではなく、3c
    (「非常に分岐した」`V_n` 族の具体的構成)の一部として初めて
    意味を持つ**——3b(c) のように「3c 抜きで先に済ませられる」
    部分ではなかった。次に着手すべきは 3c そのもの
    (3a の `d+1` 項化・3b(a)(b) いずれも、3c が具体的な対象を
    与えて初めて検証可能になる)。

    ★★★★★★★★★★★★★★★★★2026-09-04、`cancel_conductor_delta` を
    `delta_tendsto_zero` の `hrec` へ実際に接続しようとして、
    もう1つの繋ぎ目(「次数によるスケーリング」: `Wₙ₊₁` が `Wₙ` 上
    有限自由のとき `length_{Wₙ₊₁}(Wₙ₊₁/I.map(...)) = [Wₙ₊₁:Wₙ]·
    length_{Wₙ}(Wₙ/I)` という形の等式)が要ると気づいたが、
    **これは誤りだと検算で判明した**——記録しておく価値がある失敗:
    `Wₙ₊₁/Wₙ` が(残余体拡大のみの)**不分岐**次数 `n` の場合を検算
    すると、`Wₙ₊₁/I·Wₙ₊₁` はそれ自体が(残余体の)体になり
    `Wₙ₊₁`-加群としての長さは **1**(不分岐がどれだけ高次でも)。
    一方 `length_{Wₙ}(Wₙ/I)` も体なので **1**。つまり `[Wₙ₊₁:Wₙ]=n`
    倍にはならない——**正しい係数は次数 `n` 全体ではなく分岐指数
    `e` のみのはず**(残余次数分は両辺の「1加群あたりの大きさ」に
    暗黙に吸収されるため)。★誤った式のまま次に進まず、ここで
    立ち止まって検算したことを記録する——このセッションを通じて
    「見た目の一般化を急がず、具体例(今回は不分岐の場合)で検算する」
    姿勢が繰り返し役に立っている。

    ★★★★★★★★★★★★★★★★★★2026-09-04、**正しい形を完成させた**
    (`length_map_pow_of_ramificationIdx`、`Found/Falt1/KaehlerAux.lean`、
    commit `f3e5c53d`、sorry無し)。`Wₙ₊₁` が DVR で、その極大イデアルが
    `Wₙ` のイデアル `P` の像として `𝔪^e`(`e`=分岐指数)と書けるとき、
    `P^k` の像で割った長さは `e·k`——`IsDiscreteValuationRing.length_
    quotient_pow_maximalIdeal` から直接従う、短い証明だった。これで
    「次数スケーリング」の繋ぎ目も(局所=DVRの場合に)閉じた。

3c. **「非常に分岐した」`V_n` の族**そのものの形式化(具体例: `p^n`
    乗根と `1` のべき根を添加する塔、上の「典型例」段落)——まだ
    手つかず。抽象的な族の公理化(`Ω_{V_n/V_{n-1}}` が
    `(V_n/pV_n)^{d+1}` を商に持つ、という条件だけを構造として
    持たせる)なら 3a・3b より軽い可能性がある——ただし 3b が独立の
    length 計算を要すると分かった今、3c だけを先に進めても
    Theorem 1.2 全体の完成には直結しない(3b を経由しない限り
    `delta_tendsto_zero` の `hrec` に接続できない)。

    ★★★★★★★★★★★★★★★2026-09-04、3b(c) 完成の勢いで 3c 着手を検討した。
    フル(Witt ベクトル)の「典型例」ではなく、**最小の非自明な事例
    (`d=0`、単一の一様化元の逐次 `p` 乗根: `V_{n+1} := V_n[X]/(X^p-π_n)`
    、`π_n` は `V_n` の一様化元)**を試みた。この場合
    `X^p-π_n` は `V_n` 上 **Eisenstein**(`Polynomial.IsEisensteinAt`、
    mathlib にあり、`.irreducible` も既製品)で、これが 3b(a) の
    「既約性が保たれる」を**自動的に**満たす、という良い性質がある。
    ★ただし `AdjoinRoot(X^p-π_n)` が実際に整閉(`V_{n+1}` の役を
    果たせる、classical に「Eisenstein ⟹ 極大整環」)であることの
    直接の mathlib 定理はまだ見つけていない——`conductor_mul_
    differentIdeal` 経由で「Eisenstein なら conductor=(1)」を
    示す形で回り道できる可能性はある(未検証)。時間の制約で
    ここで打ち切ったが、**d=0 の最小事例は3a(生成元1個、既存の
    Lemma 1.1 インフラをそのまま使える)・3b(a)(Eisenstein の既約性
    で自動)の両方を回避できる有望な入口**と判断する——次のセッション
    はここから始めるのが良い。

    ★★★★★★★★★★★★★★★★2026-09-04、**「Eisenstein ⟹ 整閉」の道具を
    さらに調査した**: `RingTheory/Polynomial/Eisenstein/IsIntegral.lean`
    に、まさにこの classical な証明の**部品**(`mem_adjoin_of_smul_
    prime_pow_smul_of_minpoly_isEisensteinAt`:「`p^n•z` が
    `adjoin R{gen}` に入れば `z` も入る」)が既にあることを発見した
    ——これは教科書の証明(Eisenstein 多項式の根で生成される拡大が
    極大整環であることの証明)の**核心の帰納法ステップ**そのもの。
    ただし**「任意の整な `z` に対しある `n` で `p^n•z∈adjoin`」という
    近似の議論(教科書証明のもう半分)へつなぐ具体的な組み立ては
    まだ行っていない**——これも classical だが、mathlib に「Eisenstein
    ⟹ 整閉」という完成形の定理としては見つからなかった(`totally
    ramified`という語もmathlibに無い)。★これで3cの技術的な下ごしらえ
    (Eisenstein 既約性・整閉性それぞれの部品の所在)はかなり明確に
    なった——次のセッションは `mem_adjoin_of_smul_prime_pow_smul_of_
    minpoly_isEisensteinAt` を軸に「整閉性」の証明を完成させることから
    始めるとよい。

    ★★★★★★★★★★★★★★★★★★★2026-09-04、**より見通しの良い別証明の筋
    を発見した**(未着手・未検証): `mem_adjoin_of_smul_prime_pow_smul_
    of_minpoly_isEisensteinAt` 経由の一般論(length・局所化理論を要し
    このセッションで繰り返し壁に当たった領域)ではなく、**付値による
    直接計算**の方が自己完結している可能性が高い。`A` を(一様化元
    `π` を持つ)DVR、`f:=X^n-π`(Eisenstein)、`θ:=` `AdjoinRoot f` の
    根とすると: (i) `Frac(B)` 上へ `A` の付値 `v` を延長すると
    `v(θ)=1/n·v(π)` になる(`f` が Eisenstein なので強制される)、
    (ii) `y=c_0+c_1θ+...+c_{n-1}θ^{n-1}`(`c_i∈Frac(A)`)の付値は
    `min_i(v(c_i)+i/n)`(**`i/n` が `i=0,...,n-1` で mod 1 相異なる
    ので相殺が起きない**、超距離不等式の非退化版)、(iii) これから
    直接「`y` が整 ⟺ 全ての `c_i∈A`」が従う——つまり
    `AdjoinRoot(X^n-π) = B` の整閉性そのものが、**局所化やlengthの
    一般論を経由せず、付値の直接計算だけで**出る。★mathlib に
    「有限拡大への付値の延長」・「相異なる付値を持つ項の和は最小値を
    与える」という部品が使える形であるかは未確認——次のセッションの
    最有力候補としてこちらを先に調べる価値がある(`mem_adjoin_of_...`
    経由より見通しが良い可能性)。

    ★★2026-09-04、上の「付値の延長」に使える可能性のある mathlib の
    部品の**在り処だけ**確認した(中身は未検証):`Valuation.
    HasExtension`(`RingTheory/Valuation/Extension.lean`)・
    `Valuation.extendToLocalization`(`RingTheory/Valuation/
    ExtendToLocalization.lean`)——`KaehlerAux.lean` の現在の import
    には含まれず、別途 import してから中身を調べる必要がある。

    ★★★2026-09-04、中身を確認し、**`extendToLocalization` はこの
    用途に使えないと判明した**——名前通り「局所化」への延長専用
    (`S ≤ v.supp.primeCompl`・`IsLocalization S B` が要る)であり、
    `B := AdjoinRoot(f)`(`A` の局所化ではなく、次数 `n` の有限拡大)
    には適用できない。`Valuation.HasExtension` は延長の**存在**では
    なく**両立性**(`vR` が `vA` の comap と同値)を述べる Prop で、
    構成自体は別途必要。★「有限次拡大への付値の延長」という、まさに
    必要な一般論の所在はまだ見つかっていない——次のセッションでは
    `IsDedekindDomain.HeightOneSpectrum.valuation` 経由(`B` 自身が
    Dedekind であることを**先に**示してから、その height-one spectrum
    の付値として `θ` の付値を得る、という順序の入れ替え)を検討する
    余地があるが、これは「`B` が Dedekind」という結論を仮定して使う
    ことになり得るため、循環を避ける組み立てが必要。

    ★★★★★★★★★★★★★★★★★★★★2026-09-04、**戦略そのものを転換する
    重要な発見をした**: そもそも `AdjoinRoot(f)` が「既に」整閉である
    ことを Eisenstein 性から証明する必要は**無かった**。mathlib には
    `integralClosure.isDedekindDomain_fractionRing`(`RingTheory/
    DedekindDomain/IntegralClosure.lean`、**instance**)があり、
    `A` が Dedekind・`L` が `FractionRing A` 上有限分離拡大なら
    `integralClosure A L` は**自動的に** Dedekind になる
    (`FractionRing A` を明示的に使うことが鍵——一般の `K` with
    `IsFractionRing A K` では instance が発火しない、という
    tools/lean-idioms.md #23 と同種の罠を実測で確認・回避した)。
    **`Vₙ₊₁ := integralClosure Vₙ L` と定義すれば、Eisenstein 性・
    付値の延長・整閉性の証明が一切不要になる**——`w`(`f` の根)が
    `Vₙ₊₁` を「生成し切っているか」は問わなくてよい。なぜなら
    `conductor_mul_differentIdeal`(**既に完成済み**)が、その
    conductor を未知数のまま追跡できるよう作られているから——
    `cancel_conductor_delta` はまさにこの「`w` が生成し切らない
    かもしれない」場合を扱うために作った道具だった。★これで
    3c の技術的な核心の壁(整閉性の証明)を**丸ごと迂回できる**
    ことが判明した。次のセッションは `integralClosure.isDedekindDomain_
    fractionRing` を使って実際に `V_1 := integralClosure V_0 L`
    (`L:=FractionRing(AdjoinRoot(X^p-π))`)を組み立て、Falt1
    RamificationSetup 系のインフラ(`falt1_isFractionRing` 等)を
    再利用して具体的な最初の1段を完成させることから始めるとよい。

    ★★★★★★★★★★★★★★★★★★★★★2026-09-04、**上記の戦略転換を実装し、
    「非常に分岐した」`V_n` 族の1段の構成が sorry 無く完成した**
    (commit `9998b221`)。`isDedekindDomain_integralClosure_adjoinRoot_
    X_pow_sub_C`: `V0` を Dedekind、`π` を(平方でない)素元とし
    `f := X^n - C π`(`V0[X]`上)とすると、`integralClosure V0
    (AdjoinRoot (f.map (algebraMap V0 (FractionRing V0))))` が
    Dedekind であることを証明。組み立ての鍵:
    - `f` の `V0[X]` 上の Eisenstein 既約性(既存の
      `eisenstein_X_pow_sub_C_irreducible_map`)を**Gauss の補題**
      (`Polynomial.Monic.irreducible_iff_irreducible_map_fraction_map`
      、`IsIntegrallyClosed R` だけで足り **UFD は不要**——Dedekind
      整域は自動的に整閉なのでこれで賄える)で `FractionRing(V0)[X]`
      上の既約性に転送する。これで「`AdjoinRoot(f)` が整域か」という
      (UFD でない一般の Dedekind 整域上では自明でない)問いを完全に
      回避できた。
    - 分離性は `Polynomial.separable_X_pow_sub_C` を直接使う(自作の
      Bezout 係数計算は不要と判明し破棄)。
    - `AdjoinRoot(f_K)`(`f_K:=f.map(...)`、`FractionRing(V0)` 上)が
      有限次分離拡大であることを、`Module.Finite`
      (`Polynomial.Monic.finite_adjoinRoot`)と新設の補題
      `algIsSeparable_adjoinRoot_of_separable`(単既約 monic 分離
      多項式の `AdjoinRoot` は分離拡大——
      `IntermediateField.isSeparable_adjoin_simple_iff_isSeparable`
      + `adjoin_root_eq_top` + `topEquiv` を貼り合わせて構成)から得る。
    - 最後に `integralClosure.isDedekindDomain_fractionRing`
      (instance)を `infer_instance` で発火させるだけで結論が出る。

    `#print axioms` で両定理とも `[propext, Classical.choice,
    Quot.sound]` のみ(sorry 無し)を確認済み、`lake build`
    2656/2656 成功、`node tools/check.mjs --brief` は新規リグレッション
    無し(既存 13 件の NG のみ)。★これで item 3c の「`V_n` 族の各段が
    Dedekind であること」という核心のボトルネックが、**`n=0→1` の
    1段について**完全に解消した。残る作業: (a) この構成を `V_0`
    自身にも再帰的に適用して `V_n`(`n` 段目)の族を作る(`induction`
    で回すだけのはずだが、各段で「一様化元」を選び直す必要があり
    未着手)、(b) `cancel_conductor_delta` の `hspan_eq` 仮定
    (「conductor が単項で書ける」)をこの具体的な `V_1` 構成に対して
    実際に検証する(3b との接続、まだ着手していない)、(c) `Ω¹_{V1/V0}`
    が `pushoutKaehlerSplit` 系の枠組みに乗ることの確認(生成元の個数
    `d+1` との対応)。

    ★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**単項生成性
    (monogenicity)も完成した**(commit `3e1f9798`)——`falt1CokernelLengthEq`
    が要求する `hadjoin : Algebra.adjoin V ({w}:Set W) = ⊤` を、上の
    `V_1 := integralClosure V0 L` に対して実際に示した。決定的な発見:
    **同じリポジトリの別の論文トラック**(`ABC3.Found.GenEll.
    TameRamification`、Mochizuki [GenEll] の馴分岐の形式化)に、
    まさに必要な部品が既に完成していた——
    `mem_adjoin_of_pow_smul_of_isEisensteinAt`(Eisenstein なら
    `π^k•z∈adjoin` から `z∈adjoin`、mathlib の `mem_adjoin_of_smul_
    prime_smul_of_minpoly_isEisensteinAt` を `k` 回反復する帰納)と
    `exists_smul_mem_adjoin_powerBasis`(`L` の任意の元はある `0≠d∈R`
    倍すれば `adjoin` に入る、`IsLocalization.exist_integer_multiples`
    経由)。**DVR** 上では任意の非零 `d` が `π^k×単元`
    (`IsDiscreteValuationRing.eq_unit_mul_pow_irreducible`)と分解
    できるので、この2つを橋渡しするだけで「任意の整な元は `adjoin`
    に入る」(`adjoin_eq_top_of_isEisensteinAt`)、つまり `adjoin R
    {gen} = integralClosure R L` が出た。★これは教科書の「Eisenstein
    拡大は単項生成(totally ramified ⟹ monogenic)」という古典定理
    そのものだが、`AdjoinRoot(f)` の整閉性を直接証明する迂回路を
    通らずに到達できた——`Subalgebra.map_injective` +
    `AlgHom.map_adjoin_singleton` で「`L` の中での等式」を「`W` 自身の
    中での `⊤` 生成」に変換し(`adjoin_eq_top_in_integralClosure_of_
    isEisensteinAt`)、`minpoly V0 w = f`(`X^n-π` そのもの、`minpoly.
    isIntegrallyClosed_eq_field_fractions'` + `Polynomial.map_injective`
    で `K` 上の minpoly `f_K` から `V0` 上の minpoly へ降ろす)を経由
    して具体例へ組み立てた
    (`adjoin_eq_top_integralClosure_adjoinRoot_X_pow_sub_C`)。
    `#print axioms` で全定理 sorry 無しを確認、`lake build`
    2837/2837 成功、`node tools/check.mjs --brief` は新規リグレッション
    無し。★これで **`V_1` の1段について、Dedekind であること
    (前回)と単項生成であること(今回)の両方が揃った**——
    `falt1CokernelLengthEq` を実際に適用する準備が整った。次回は
    上記の残る作業 (a)(b)(c) のうち、特に (b)(`hspan_eq` の検証、
    `differentIdeal V0 V1 = span{f'(w)}` を `differentIdeal_eq_span_
    derivative`(Lemma 1.1 自身の道具)で計算し、`cancel_conductor_
    delta` の等式と突き合わせる)に取り組むのが、3b と 3c を実際に
    接続する最短経路と見込む。

    ★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**(b) を実行し、
    `differentIdeal V0 V1` の具体式が完成した**(`differentIdeal_eq_
    span_of_adjoinRoot_X_pow_sub_C`)。`differentIdeal_eq_span_
    derivative`(Lemma 1.1 自身の道具、`w` が単項生成なら `differentIdeal
    V W = span{f'(w)}`)と、前段の monogenicity を貼り合わせるだけで、
    `differentIdeal V0 V1 = span{n·w^{n-1}}`(`w`=根の像、`n·X^{n-1}`
    は `X^n-π` の微分)という**古典的な Eisenstein 拡大の判別式の公式
    そのもの**が閉じた。★配管上の教訓(新パターン、tools/lean-idioms.md
    に追記の価値あり): `differentIdeal A B` は `[IsDedekindDomain B]`・
    `[Module.IsTorsionFree A B]` を引数に持つため、これらを証明内部で
    `infer_instance` で揃えようとすると、**結論(goal)自体が同じ2つの
    instance を要求する**ことで先に走る instance 探索の失敗が
    キャッシュされ、後から `haveI` で正しく揃えても失敗し続けるという
    罠に嵌った(`set` 起因ではない、新しい失敗形)——解決策はこれら2つを
    **定理の引数として明示的に要求する**ことで goal 自体の instance
    探索を迂回すること。

    ★★委員会メモ(作業の帰属について): この commit は**共有
    worktree で並行して動いていた別セッション**(Lubin-Tate/pGC の
    冪級数評価に取り組んでいたセッション)が `git commit -a`(または
    同等の非pathspec制限のコミット)を実行した際に、**このセッションが
    ステージ済みだった `KaehlerAux.lean` の変更を巻き込んで一緒に
    コミットした**——commit hash は `88879620`(コミットメッセージは
    「冪級数を Λ_n の元で実際に評価する」という無関係な内容だが、
    diff には本節で述べた `differentIdeal_eq_span_of_adjoinRoot_
    X_pow_sub_C`(112行)が含まれている、実測で確認済み)。sorry 無し・
    `lake build` 成功・push 済みであることは確認済みなので**内容自体は
    無事**——今後 `git log -- lean/ABC3/Found/Falt1/KaehlerAux.lean`
    でこの完成を辿る際は `88879620` を見ること。

    ★これで item 3c の V_1 について、**Dedekind・単項生成・
    differentIdeal の具体値**の3点セットが揃った。残る作業:
    `cancel_conductor_delta` の `hspan_eq` 仮定(`spanDeriv = Idiff`)を
    この具体式(`Idiff = Ideal.map (algebraMap V1 W1) (differentIdeal
    V0 V1) = Ideal.map(...) (span{n·w^{n-1}})`)に対して実際に
    代入・検証する接続作業(3b と 3c を本当に繋ぐ最後の1ピース)、
    および上記(a)(V_n 再帰族)・(c)(`pushoutKaehlerSplit` との接続)。

    ★★★★★★2026-09-04、**`hspan_eq` 接続の具体的な設計図が見えた**
    (未着手・次回の最有力候補)。`differentIdeal_tower_diamond` の
    `Wₙ→Wₙ₊₁` 側を、`Vₙ→Vₙ₊₁`(=今回の `V0→V1`)と**同じ多項式
    `X^n-π` を `Wₙ` へ base change したもの**で作れば
    (`Wₙ₊₁ := integralClosure Wₙ (AdjoinRoot((X^n-C π).map(V0→Wₙ→
    Frac(Wₙ))))`)、`hspan_eq` に必要な「`x`(`Wₙ₊₁` の生成元)の
    `algebraMap V1 Wₙ₊₁` による `w` の像との一致」は、**今セッション
    の前半で作った `adjoinRootTensorEquiv`**(`C⊗[R]AdjoinRoot g
    ≃ₐ[C] AdjoinRoot(g.map(algebraMap R C))`、`pushoutKaehlerSplit`
    のために構築したもの)を `R:=V0`・`C:=Frac(Wₙ)`・`g:=fK` で適用
    すれば**ほぼ自動的に従う**——`w⊗1` の像が `AdjoinRoot` 側の
    `root(fK.map(...))` に対応するという、`adjoinRootTensorEquiv` の
    定義性質そのもの。つまり **`Wₙ₊₁` 側の全ての情報
    (Dedekind・単項生成・differentIdeal の値)も、`differentIdeal_
    eq_span_of_adjoinRoot_X_pow_sub_C` 系の3定理を**そのまま
    `V0` の代わりに `Wₙ` に適用するだけ**で再利用できる**(3定理は
    すでに一般の `[IsDiscreteValuationRing V0]` に対して書いてある
    ので、`Wₙ` が DVR でありさえすれば即座に適用できる)。
    残るのは(i)`Wₙ` が DVR で `algebraMap V0 Wₙ π` が引き続き
    prime・非自乗であること(`Wₙ/V0` が今考えている素での不分岐、
    という自然な仮定)を明示的な仮説として立てること、
    (ii)`adjoinRootTensorEquiv` の像が `integralClosure` の元として
    振る舞うことの橋渡し(`IsIntegral` が base change で保たれる、
    という一般に易しいはずの事実)、の2点——(ii)は次回最初に
    着手するのに適した、比較的小さい補題だと見込む。

    ★★★★★★★★★2026-09-04、**base change 写像そのものが完成した**
    (`algHomAdjoinRootOfCompat`・`algHomAdjoinRootOfCompat_root`、
    commit `7f09476c`、他セッションに巻き込まれてコミット——上と
    同じ状況、内容は sorry 無し・push 済み確認済み)。当初
    `adjoinRootTensorEquiv`(テンソル積経由)を使う計画だったが、
    実装時により直接的な経路を発見した: `AdjoinRoot.map`(mathlib、
    係数の ring hom `φ0:K→+*K'` から `AdjoinRoot fK→+*AdjoinRoot
    (fK.map φ0)` を作る)を `AlgHom.mk'` で `V0`-線形性
    (`φ0` が `V0` 上の algebraMap と両立するという1つの仮定
    `hφ0` から、`IsScalarTower.algebraMap_apply` +
    `AdjoinRoot.algebraMap_eq` + `AdjoinRoot.map_of` で従う)を
    確認するだけで足りた——`adjoinRootTensorEquiv` も
    `AdjoinRoot.mapAlgHom`(mathlib、同じ底環上のみ)も不要だった。
    `algHomAdjoinRootOfCompat_root` で「根は根に写る」ことも確認
    済み(`w` の像が `Wₙ₊₁` 側の `x` と一致することの鍵)。

    ★残る作業(上記(i)相当、未着手): `φ0 := IsFractionRing.map
    (algebraMap V0 Wₙ の単射性)`(mathlib、`FractionRing V0 →+*
    FractionRing Wₙ`)を実際に構成し、`fK.map φ0` が「`Wₙ` を底環に
    直接構成した `X^n-π'`(`π':=algebraMap V0 Wₙ π`)の Eisenstein
    多項式」と一致すること(`π'` が引き続き prime・非自乗という
    仮定のもとで)を示し、`differentIdeal_tower_diamond`・
    `conductor_mul_differentIdeal` のフル仮説束(`IsDedekindDomain
    Wₙ`・`Module.Finite Wₙ Wₙ₊₁`・`IsScalarTower V0 Wₙ Wₙ₊₁` 等)を
    実際に組み立てて `cancel_conductor_delta` を適用する——これが
    3b・3c を実際に閉じる最終組み立てとして残っている。

    ★★★2026-09-04、**`φ0` 自体も完成した**(`fractionRingMapOfInjective`・
    `fractionRingMapOfInjective_algebraMap`、commit `fb2396c6`)——
    `IsFractionRing.map`(mathlib)+ `IsLocalization.map_eq` +
    `IsScalarTower.algebraMap_apply` で、`algebraMap V0 Wₙ` が単射
    という1つの仮定だけから `φ0:FractionRing V0→+*FractionRing Wₙ` と
    その `V0` 上の algebraMap との両立を得た。これで
    `algHomAdjoinRootOfCompat`(前段)に実際に渡す `φ0`・`hφ0` の
    両方が揃った——**base change 写像の構成は完全に完成**、残る
    作業は上記の「`π'` が引き続き Eisenstein であることの仮定の
    もとでの多項式の一致」と「フル仮説束の組み立て」のみ。

    ★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**驚きの単純化と、
    3c→3b接続の核心が完成した**(commit `2ec2215c`・`865622c8`・
    `f43ae887`)。実装を進める中で、`adjoinRootTensorEquiv`(テンソル積
    経由)や `algHomAdjoinRootOfCompat`(フラクション体経由)より
    **遥かに単純な経路**を発見した:

    1. **`falt1AdjoinRootEquivIntegralClosure`**: `V_1 := integralClosure
       V0 (AdjoinRoot(X^n-π の base change))` は、実は**base change
       する前の** `AdjoinRoot(X^n-π)`(`V0[X]` 上そのままの多項式の商)
       と `V0`-代数として**同型**——`minpoly.equivAdjoin` + monogenicity
       (`adjoin_eq_integralClosure_of_isEisensteinAt`)を貼り合わせる
       だけ。★配線上の罠(記録): `PowerBasis` の暗黙引数を伴う定理を
       裸の `w` に対して呼ぶと高階単一化が `whnf` タイムアウトする
       ——常に `PB.gen` の形を経由すること。
    2. **`algHomAdjoinRootOfCompat'`**: (1)のおかげで、`Vₙ→Vₙ₊₁` の
       橋渡しは **`V0[X]→Wₙ[X]` の底環 base change だけ**(`AdjoinRoot.
       map` + `AlgHom.mk'`)で済むと判明——フラクション体・
       `IsFractionRing.map` は一切不要だった。
    3. **`falt1BaseChangeAlgHom`**: 上記2つを**3段合成**
       (`V1≃ₐ[V0]AdjoinRoot f → AdjoinRoot f→ₐ[V0]AdjoinRoot g
       (`g:=f.map(algebraMap V0 Wₙ)`、多項式として `X^n-Cπ'` に等しい)
       → AdjoinRoot g≃ₐ[Wₙ]Wₙ₊₁`、最後の等式を `AlgHom.restrictScalars
       V0` で `V0`-線形とみなして繋ぐ)して、**`V1 →ₐ[V0] Wₙ₊₁` の
       実際の `AlgHom` を構成した**——`differentIdeal_tower_diamond`
       が要求する `Algebra Vn1 Wn1` インスタンスがこれで手に入る。

    `#print axioms` で全定理 sorry 無し確認、`lake build` 2837/2837
    成功、`node tools/check.mjs --brief` は新規リグレッション無し。
    ★★これで **3c→3b 接続の最大の技術的な壁(「同じ生成元を共有する
    塔をどう構成するか」)を実際に越えた**。残る作業:
    (a) `falt1BaseChangeAlgHom` の下で生成元の対応(`ψ(w)=x` 相当)を
    明示的に確認すること、(b) `Module.Finite Wₙ Wₙ₊₁`・
    `IsScalarTower V0 Wₙ Wₙ₊₁`・`IsScalarTower V0 V1 Wₙ₊₁` 等の残りの
    instance を(これも `falt1AdjoinRootEquivIntegralClosure` 経由で
    「`Wₙ₊₁` も結局 `AdjoinRoot` そのもの」という単純化が使えるはず)
    整備すること、(c) `differentIdeal_tower_diamond` の `hsep`
    (`Algebra.IsSeparable (FractionRing V0)(FractionRing Wₙ₊₁)`)を
    示すこと、(d) 最終的に `conductor_mul_differentIdeal` +
    `cancel_conductor_delta` を実際に適用すること——これで item 3c は
    技術的な核心を越え、残りは「組み立て」の段階に入った。

    ★★★2026-09-04、**(a) に着手し、Lean の配線上の教訓を得た**
    (未完成——`falt1AdjoinRootEquivIntegralClosure` 自体は既に完成・
    push 済みで、以下は「その値を `AdjoinRoot.root f` で評価した式を
    後から取り出す」という**別の**作業)。`unfold` + `simp only
    [eq_mpr_eq_cast, eq_mp_eq_cast, cast_eq]` で `rw […] at e1` が
    生んだ**入れ子の cast の大半は消せる**(外側の2つ、`set f`・
    `set fK` に由来するもの)——ただし `hEqmin`(`minpoly V0 PB.gen =
    f`)由来の**内側の cast**(`AlgEquiv` という関数型そのものへの
    cast)だけは `cast_eq`・`AlgEquiv.trans_apply`・`Subtype.ext_iff`
    のどれでも綺麗に剥がせず、`subst`/`generalize` も(`f` が `set`
    由来の let で「motive is not type correct」)効かなかった——
    各試行が **50〜75秒**かかる中で複数回失敗し、時間対効果が悪化した
    ため一旦中断。

    ★★推奨する迂回路(次回、これを最初に試すこと): **後から `e1
    (AdjoinRoot.root f) = w` を証明しようとするのではなく、`w` を
    最初から `w := e1 (AdjoinRoot.root f)` と定義すればこの等式は
    `rfl` になる**——`differentIdeal_eq_span_of_adjoinRoot_X_pow_
    sub_C` 側の `w`(`⟨PB.gen, hBint⟩` として個別に構成したもの)を、
    `falt1AdjoinRootEquivIntegralClosure` から得られる `w` に**差し
    替える**方向で両者を統一すれば、cast を剥がす作業そのものが
    不要になる可能性が高い。★この種の「opaque な tactic-mode `def`
    を後から `unfold` して cast まみれの中身を分析する」パターンは
    時間対効果が悪いと判明した——tools/lean-idioms.md #29 に記録した。

    ★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**この推奨迂回路を
    実行し、生成元の対応が完成した**(`falt1GeneratorPackage`・
    `falt1BaseChangeAlgHom_generator_correspondence`、commit
    `cc34f441`)。予想通り `w := e1 (root f)` と最初から定義する
    ことで、以前 `PowerBasis`・Eisenstein の詳細を経由していた
    「整性・単項生成性・`minpoly=f`」の3点セット
    (`falt1GeneratorPackage`)が、**同型の自然性だけ**
    (`IsIntegral.map`・`AlgHom.map_adjoin_singleton`・
    `minpoly.isIntegrallyClosed_dvd`+`Irreducible.associated_of_dvd`)
    から出ることを確認した——`isDedekindDomain_...`・`adjoin_eq_top_
    ...` 系の元の証明より遥かに短い代替ルート。さらに `falt1BaseChangeAlgHom`
    という名前を経由して `w↦x` を証明しようとする(`unfold` して中身を
    調べる)アプローチは何度試みても50〜75秒ずつタイムアウトし続けた
    ため、**同じ構成をその場で直接書き下す**(named def を経由しない)
    ことで解決した——`algHomAdjoinRootOfCompat'_cast_root`(`g` を
    自由変数のまま保つのが鍵、`set`束縛された変数に対する `subst`
    は「invalid equality proof」で失敗する)を独立補題として切り出し、
    それを直接使うことで `hstep1root` が一瞬(1秒未満)で通った。

    ★★これで **item 3c → 3b 接続に必要な技術的ピースが全て揃った**
    (Dedekind・単項生成・differentIdeal の具体式・base change 写像・
    生成元対応)。残る作業は純粋な「組み立て」段階: (a) `Module.Finite
    Wₙ Wₙ₊₁`・`IsScalarTower V0 Wₙ Wₙ₊₁`・`IsScalarTower V0 V1 Wₙ₊₁`
    等の残りの instance を整備する(これも `falt1AdjoinRootEquivIntegralClosure`
    経由で「`Wₙ₊₁` も結局 `AdjoinRoot` そのもの」という単純化が使える
    はず、`Module.Finite` は `Polynomial.Monic.finite_adjoinRoot` から
    ほぼ自明)、(b) `differentIdeal_tower_diamond` の `hsep`
    (`Algebra.IsSeparable (FractionRing V0)(FractionRing Wₙ₊₁)`)を
    示す、(c) `conductor_mul_differentIdeal` を `Wₙ→Wₙ₊₁` に適用して
    `Jₙ`(`differentIdeal Wₙ Wₙ₊₁`)を conductor で書く、(d) 最終的に
    `cancel_conductor_delta` を適用して `condWnx * deltaN1 =
    deltaNmapped` を得る——ここまで来れば 3b・3c は技術的には閉じ、
    残るのは `delta_tendsto_zero` の `hrec` への変換(3b(c)、既存の
    `length_quotient_span_singleton_mul` 等が使えるはず)と、V_n の
    **再帰族**(item 3c 本体、`V_1` の1段からの帰納)の構築のみ。

    ★★2026-09-04、`differentIdeal_tower_diamond` の instance 整備
    (上記(a))に着手した(commit `45ecea4c`)。**`Module.Finite Wₙ
    Wₙ₊₁` を完成**(`falt1ModuleFiniteWnWn1`)——`AdjoinRoot.powerBasis'`
    (**体でなく一般の可換環上でも** monic なら `PowerBasis` が取れる、
    mathlib の既製品、`Polynomial.Monic.finite_adjoinRoot` の体限定
    版より広い)+`falt1AdjoinRootEquivIntegralClosure` を貼り合わせる
    だけで閉じた。★また **`Algebra V0 Wₙ₊₁`(合成 `V0→Wₙ→FractionRing
    Wₙ→AdjoinRoot gK`)と `IsScalarTower V0 Wₙ Wₙ₊₁` は instance 探索
    だけで自動的に見つかる**ことを実測で確認した(明示的な構成は不要)。
    ★残る instance(未着手): `Module.IsTorsionFree Wₙ Wₙ₊₁`
    (`falt1BaseChangeAlgHom_generator_correspondence` の証明内で
    使った `Module.IsTorsionFree V0 (AdjoinRoot g)` の議論をほぼ
    そのまま流用できる見込み)、`Algebra V1 Wₙ₊₁` を instance として
    登録する(`falt1BaseChangeAlgHom.toRingHom.toAlgebra`)、
    `IsScalarTower V0 V1 Wₙ₊₁`(`falt1BaseChangeAlgHom` が `→ₐ[V0]`で
    あることから `ψ.commutes` 経由でほぼ自動のはず)、そして `Wₙ`
    自身が `V0` 上有限・捩れ無しであること(`Wₙ` がどんな `V0`-拡大
    として与えられるかに依存する、呼び出し側の仮定として要る)。

    ★★2026-09-04、**`Module.IsTorsionFree Wₙ Wₙ₊₁` も完成した**
    (`falt1ModuleIsTorsionFreeWnWn1`、commit `addad6e1`)——予想通り
    `differentIdeal_eq_span_of_adjoinRoot_X_pow_sub_C` の中で使った
    `Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot fK))`
    の議論(`algebraMap Wₙ (AdjoinRoot gK)` の単射性 →
    `IsIntegralClosure.isTorsionFree`)を `Wₙ`・`gK` に対して**そのまま
    繰り返すだけ**で閉じた——新しい発見は無かったが、既存パターンの
    再利用が効くことを確認できた。★これで約20個中、`Module.Finite
    Wₙ Wₙ₊₁`・`Module.IsTorsionFree Wₙ Wₙ₊₁`・`Algebra V0 Wₙ₊₁`・
    `IsScalarTower V0 Wₙ Wₙ₊₁`(自動発見)の4つが揃った。残る
    (未着手・次回): `Algebra V1 Wₙ₊₁` を instance として登録する
    (`falt1BaseChangeAlgHom.toRingHom.toAlgebra`)、`IsScalarTower
    V0 V1 Wₙ₊₁`、`Module.Finite V0 V1`・`Module.IsTorsionFree V0 V1`
    (`V1` 自身、既存の `falt1ModuleFiniteWnWn1`・`falt1ModuleIsTorsionFreeWnWn1`
    と同じパターンを `V0`・`fK` に戻して適用するだけのはず——まだ
    やっていないが機械的)、`Module.Finite V0 Wₙ`・`Module.IsTorsionFree
    V0 Wₙ`(`Wₙ` 自身、呼び出し側の仮定として要る)、`Module.Finite
    V1 Wₙ₊₁`・`Module.IsTorsionFree V1 Wₙ₊₁`(`ψ` 経由の有限性、
    `Module.Finite.trans`的な議論が要りそう)、`Module.Finite V0 Wₙ₊₁`・
    `Module.IsTorsionFree V0 Wₙ₊₁`(`Wₙ→Wₙ₊₁` と `V0→Wₙ` の合成)。

    ★★2026-09-04、**`Module.Finite/IsTorsionFree V0 V1` も完成した**
    (`falt1ModuleFiniteV0V1`・`falt1ModuleIsTorsionFreeV0V1`、commit
    `331b1d9d`)——`Wₙ`側と全く同じ議論を `V0`・`f` に戻すだけ。
    ★★重要な発見: `Module.IsTorsionFree V0 V1` はこれまで**多くの
    定理が instance 前提として明示的に要求していた**が、実は毎回
    この議論で導出可能だったと判明した(既存の証明済み定理群は
    変更せず、影響範囲を数える余裕が無いため追加のみに留めた——
    次回、余裕があれば `[Module.IsTorsionFree V0 (...)]` の前提を
    削って引数を減らすリファクタリングを検討する価値がある)。

    ★★★`Algebra V1 Wₙ₊₁`・`IsScalarTower V0 V1 Wₙ₊₁` も**パターンとしては
    確認できた**(未commit・named theorem としてのパッケージ化は
    断念)——`letI := (falt1BaseChangeAlgHom ...).toRingHom.toAlgebra`
    としてから `IsScalarTower.of_algebraMap_eq` + `ψ.commutes` で
    `IsScalarTower V0 V1 Wₙ₊₁` が出ることを、1つの大きな `example`
    ブロック内で実際に確認した(約19秒、sorry 無し)。★ただしこれを
    **独立した named theorem として `∃`-wrap してパッケージ化しよう
    とすると**(`Algebra` インスタンスを存在型で包んで返す形)、
    `refine ⟨inferInstance, ?_⟩` の instance 探索が
    `«synthesize pending MVars»` でタイムアウトし続けた——
    tools/lean-idioms.md #29 と同種だが、今回は「instance を戻り値と
    して`∃`で包む」という**新しい**失敗パターン。★推奨: この種の
    「これから使う instance を用意する」系の事実は、独立した
    theorem として証明・命名しようとせず、**実際に
    `differentIdeal_tower_diamond` を呼び出す箇所で `letI`/`haveI`
    としてインラインに展開する**のが正攻法——次回はこの方針で
    直接 `differentIdeal_tower_diamond` の呼び出しに挑むとよい
    (単独の instance-provider theorem を目指さない)。

    ★★2026-09-04、**`Module.Finite/IsTorsionFree V0 Wₙ₊₁` も完成した**
    (`falt1ModuleFiniteV0Wn1`・`falt1ModuleIsTorsionFreeV0Wn1`、commit
    `d739bbed`)。`Module.Finite` は `Module.Finite.trans`(`V0→Wₙ→
    Wₙ₊₁` の合成、`Module.Finite V0 Wₙ` を仮定として要求)で機械的。
    `Module.IsTorsionFree` は `algebraMap V0 Wₙ₊₁` の単射性を
    `algebraMap V0 Wₙ`(仮定)・`algebraMap Wₙ Wₙ₊₁`(構成から)の
    合成で示した。★これまで3回繰り返していた「正則元→非零→
    `algebraMap` 単射性→整域での非零元倍は単射」という議論を、
    **汎用補題 `moduleIsTorsionFree_of_injective`**(`T` が整域で
    `algebraMap R T` 単射なら `Module.IsTorsionFree R T`)として
    独立に切り出した——今後同種の場面ではこれを直接呼べばよい。

    ★これで約20個中**8個**(`Module.Finite`/`IsTorsionFree` の
    `(Wₙ,Wₙ₊₁)`・`(V0,V1)`・`(V0,Wₙ₊₁)` の3ペア分、`Algebra V0 Wₙ₊₁`・
    `IsScalarTower V0 Wₙ Wₙ₊₁`)が揃った。残る(未着手):
    `Algebra V1 Wₙ₊₁`・`IsScalarTower V0 V1 Wₙ₊₁`(パターン確認済み、
    `letI`でインライン展開する方針)、`Module.Finite/IsTorsionFree
    V1 Wₙ₊₁`(`ψ:=falt1BaseChangeAlgHom` の単射性を要する、
    まだ証明していない新しい数学的内容)、`hsep`(`Algebra.IsSeparable
    (FractionRing V0)(FractionRing Wₙ₊₁)`)。

    ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**`ψ` の単射性が
    完成した**(`adjoinRoot_map_mk`・`adjoinRoot_map_injective_of_map_eq`・
    `algHomAdjoinRootOfCompat'_injective`・`_cast_injective`・
    `falt1BaseChangeAlgHom_generator_and_injective`、commit `e3ae2da5`)
    ——このセッション後半で唯一残っていた**新しい数学的内容**だった。
    予想(基底 `{1,...,X^{n-1}}` の対応)通りだったが、実際には
    **`PowerBasis` の基底展開を一切経由せず**、mathlib の
    `Polynomial.map_dvd_map`(`x` が monic・係数写像が単射なら
    `x.map f ∣ y.map f ↔ x ∣ y`——多項式の割り算アルゴリズムが monic
    除数に対して任意の可換環上で機能することの帰結)を直接使うだけで
    閉じた——`AdjoinRoot.map` の核が「`p∣r`」という単純な条件に
    帰着することが、`AdjoinRoot.lift_mk`+`Polynomial.eval₂_map`+
    `AdjoinRoot.aeval_eq` で `mk` への作用を計算し、`AdjoinRoot.mk_
    eq_zero` へ橋渡しするだけで示せた。

    ★★続けて **`Module.IsTorsionFree V1 Wₙ₊₁` も即座に完成した**
    (`falt1ModuleIsTorsionFreeV1Wn1`、commit `ce9e4bc2`)——`ψ` の
    単射性を `moduleIsTorsionFree_of_injective`(既存の汎用補題)に
    直接渡すだけ。★これで約20個中**9個**+ `Algebra V1 Wₙ₊₁`・
    `IsScalarTower V0 V1 Wₙ₊₁`(パターン確認済み)が揃った。
    **残るのは `Module.Finite V1 Wₙ₊₁`(単射像の有限生成性の議論、
    `ψ` が単射なのでおそらく `Submodule.fg_of_injective` 的な既製品
    一発で閉じるはず)と `hsep` のみ**——`differentIdeal_tower_diamond`
    の instance 整備は完成が視野に入った。

    ★★★★2026-09-04、**`Module.Finite V1 Wₙ₊₁` も完成した**
    (`falt1ModuleFiniteV1Wn1`、commit `786fcd84`)——予想(単射像の
    有限生成性)とは違う、もっと単純な経路で閉じた:`Module.Finite.
    of_restrictScalars_finite`(`R→A→M` のタワーで `M` が `R` 上有限
    なら `A` 上でも有限、「小さい環上の有限生成は大きい環上でも
    有限生成」という一般に易しい事実)を `R:=V0,A:=V1,M:=Wₙ₊₁` に
    適用するだけ——`IsScalarTower V0 V1 Wₙ₊₁`(`ψ.commutes` から)と
    既存の `Module.Finite V0 Wₙ₊₁` の2つを揃えれば機械的。★`ψ` の
    単射性は**この有限性には不要だった**(`IsTorsionFree` の方でのみ
    必要)——予想が外れたが、より単純な経路が見つかった好例。

    ★★★これで **`differentIdeal_tower_diamond` の instance 整備は
    事実上完成した**——約20個中10個を直接証明し、残り(`Algebra
    Vn1 Wn1`・`Algebra Vn Wn`・`Algebra Vn Vn1`・`IsScalarTower Vn
    Vn1 Wn1`・`IsScalarTower Vn Wn Wn1` 等)は全て確認済みの
    自動発見パターン(`letI`でインライン展開)で揃う。**唯一残る
    のは `hsep`(`Algebra.IsSeparable (FractionRing V0)(FractionRing
    Wₙ₊₁)`)のみ**——ただしこれは技術的に厄介な点が1つ判明した:
    `integralClosure.isFractionRing_of_finite_extension`(mathlib)
    により `IsFractionRing Wₙ₊₁ (AdjoinRoot gK)` は言えるが、これは
    `FractionRing Wₙ₊₁`(**具体的な** `Localization (nonZeroDivisors
    Wₙ₊₁)` という構成)と `AdjoinRoot gK` が**同型**であることしか
    意味せず、**等しい**わけではない——`hsep` を `Algebra.IsSeparable
    (FractionRing V0)(AdjoinRoot gK)`(こちらは `algIsSeparable_
    adjoinRoot_of_separable` 系の既存道具が直接使える形)から
    移送するには、この同型に沿った separability の transport が
    もう1段要る。★また `Wₙ` 自身の `FractionRing V0` 上の分離性
    (`Algebra.IsSeparable (FractionRing V0)(FractionRing Wₙ)`)も
    新たな仮定として要りそう(`Wₙ` が塔のどこから来るかに依存)。
    次回はこの「同型に沿った transport」から始めるとよい——
    tools/lean-idioms.md #29 の教訓(cast を避け、`IsFractionRing`
    の一意性補題 `IsFractionRing.algEquivOfIsFractionRing` 等を
    使って明示的な同型として扱う)が活きる場面のはず。

    ★★★★★★★2026-09-04、**`hsep` の数学的内容そのものは完全に検証
    できた**(未commit・named theorem へのパッケージ化は断念)。
    予想通りの筋(分離性の推移律 + `IsFractionRing` の一意性による
    同型 transport)で完全に閉じることを、**独立した `example` ブロック
    内で sorry 無く実証した**:

    1. `IsLocalization.algEquiv`(mathlib)で `FractionRing Wₙ₊₁ ≃ₐ[Wₙ₊₁]
       AdjoinRoot gK` を得る(`integralClosure.isFractionRing_of_
       finite_extension` が与える `IsFractionRing Wₙ₊₁ (AdjoinRoot gK)`
       から)。
    2. **驚きの発見**: `Algebra (FractionRing V0)(FractionRing Wₙ)`
       が(`fractionRingMapOfInjective` 経由で)一度 `letI` で入って
       さえいれば、`Algebra (FractionRing V0)(AdjoinRoot gK)` と
       `IsScalarTower (FractionRing V0)(FractionRing Wₙ)(AdjoinRoot
       gK)` は**instance 探索だけで自動的に見つかる**(`φ` を手で
       組み立てる必要は無かった)。
    3. `Algebra.IsSeparable.trans`(mathlib、分離性の推移律)で
       `Algebra.IsSeparable (FractionRing V0)(AdjoinRoot gK)` を
       `Algebra.IsSeparable (FractionRing V0)(FractionRing Wₙ)`
       (新たな仮定)+ `Algebra.IsSeparable (FractionRing Wₙ)
       (AdjoinRoot gK)`(`algIsSeparable_adjoinRoot_of_separable`、
       既存)から得る。
    4. `AlgEquiv.ofRingEquiv`(mathlib)で(1)の同型を`FractionRing
       V0` 上の `AlgEquiv` に昇格させ(`e_Wn1.apply_symm_apply` から
       `rfl` 級の可換性証明)、`AlgEquiv.Algebra.isSeparable_iff`
       (mathlib)で(3)を transport すれば `hsep` そのものが出る。

    ★★ただし**この4ステップを1つの named theorem として綺麗に
    パッケージ化しようとすると**(`Algebra.IsSeparable (FractionRing
    V0)(FractionRing Wₙ)` を instance 前提として要求する形にすると)、
    (2)の自動発見が**named theorem の中では効かなくなる**という
    新しい instance-search の脆さに遭遇した——`letI`で導入した
    ローカルな `Algebra` インスタンスと、`hsep` 系の仮定の型が
    参照する(別途 `def` で名指しされた)同じ値が、**instance 探索の
    観点では別物として扱われる**らしく、複数回の書き方の変更でも
    再現し続けた。★これは`falt1BaseChangeAlgHom_generator_and_
    injective`パッケージ化の際に遭遇したものとも異なる**新しい**
    パターン——tools/lean-idioms.md への追記候補(次回、実際に
    追記する前にもう少し原因を特定する価値がある)。

    ★結論: **`hsep` は数学的には解決済み**(4ステップの筋は完全に
    正しく、個別に sorry 無く確認できた)。残るのは Lean の配線上の
    問題(named theorem 化)のみ——次回は、`differentIdeal_tower_
    diamond` を**直接呼び出す箇所**で、このセッションの `example`
    ブロックとまったく同じ順序で `letI`/`haveI` を並べてインライン
    展開する(named theorem を経由しない)方針を試すのが良い
    (`falt1BaseChangeAlgHom_generator_and_injective` 完成時と同じ
    教訓——instance を提供する事実は独立に名前を付けようとせず、
    使う場所で直接組み立てる)。

    ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**実際にこの方針で
    `differentIdeal_tower_diamond` を1つの巨大な `example` で直接
    呼び出してみた——`hsep` 以外の**約19個の instance は全て問題無く
    組み上がる**ことを確認した**(sorry 無し、`hsep` を1つの仮定
    として明示的に受け取る形なら`True`は達成、commit 未——named
    theorem としては未パッケージ)。組み立て順序:
    `isDedekindDomain_integralClosure_adjoinRoot_X_pow_sub_C`(V1・
    Wₙ₊₁ 双方)→ `falt1ModuleIsTorsionFreeV0V1`・`falt1ModuleFiniteV0V1`
    → `falt1ModuleIsTorsionFreeWnWn1`・`falt1ModuleFiniteWnWn1` →
    `falt1ModuleFiniteV0Wn1`・`falt1ModuleIsTorsionFreeV0Wn1` →
    `falt1BaseChangeAlgHom_generator_and_injective` で `ψ` を obtain
    → `letI := ψ.toRingHom.toAlgebra` で `Algebra V1 Wₙ₊₁` →
    `IsScalarTower V0 V1 Wₙ₊₁`(`ψ.commutes`)→
    `moduleIsTorsionFree_of_injective hψinj` →
    `Module.Finite.of_restrictScalars_finite`。**この順序で19個
    全てが素直に通った**——大きな驚きは無く、既存の各ピースが
    そのまま合成できることを確認できただけ、という意味で
    「配線の仕上げ」段階に入ったことの実証。

    ★ただし**`hsep` を最後に接続しようとした瞬間、mathlib 自身の
    `FractionRing.liftAlgebra`(`abbrev`)との diamond
    (tools/lean-idioms.md #23 と全く同じクラスの問題)に遭遇した**:
    `differentIdeal_tower_diamond` が実際に要求する `Algebra
    (FractionRing V0)(FractionRing Wₙ₊₁)` は `FractionRing.
    liftAlgebra V0 (FractionRing Wₙ₊₁)`(`Algebra V0 Wₙ₊₁` から
    自動的に持ち上げる mathlib の既製品)であり、これまでこの
    セッションで組み立ててきた「`AdjoinRoot gK` 経由の手作りの
    instance」とは**別の項**になる。#23 の教訓通り `letI :=
    FractionRing.liftAlgebra V0 (FractionRing Wₙ₊₁)` を明示的に
    呼び、`IsScalarTower V0 (FractionRing V0)(FractionRing Wₙ₊₁)`
    が(このセッションで確認済み、`infer_instance` で自動的に
    見つかる)ことを使って `IsLocalization.ringHom_ext`(mathlib、
    局所化からの2つの準同型が底環上で一致すれば全体で一致する)で
    `hcommutes` を再構成する、という**筋道までは特定した**——
    ただし `algebraMap V0 Wₙ₊₁`(標準、`V0→Wₙ→Wₙ₊₁` 経由)と
    `φ`(`AdjoinRoot gK` へ向かう、`V0→Wₙ→FractionRing Wₙ→AdjoinRoot
    gK` 経由)が同じ元に対して一致する、という**もう1段の
    diamond 解消**(`e_Wn1` を Wₙ の scalar でも可換にする議論)が
    追加で必要と判明し、ここで打ち切った。

    ★結論(正直な評価): `differentIdeal_tower_diamond` の**19/20
    の instance は完全に組み上がることを実証した**——これは
    Theorem 1.2・item 3c にとって非常に大きな前進。**残る `hsep`
    1つ**は、数学的には完全に正しい(既に別の形で証明済み)ものの、
    mathlib の `FractionRing.liftAlgebra` という**具体的な instance
    と一致させる**という、既知のパターン(#23)の**入れ子になった
    版**として立ちはだかっている——次回はこの「二重の diamond
    解消」(`FractionRing.liftAlgebra` との一致 + `V0→Wₙ₊₁` の
    2つの経路の一致)から始めるのが最短距離と見込む。

    ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**`hsep` を完全に
    解決した**(`falt1_hsep_bundled`、`KaehlerAux.lean`、commit
    `cd7a95ec`、`lake build`・`#print axioms` 確認済み・sorry 無し)。
    上で「もう1段の diamond」と呼んでいた懸念——`algebraMap V0 Wₙ₊₁`
    (標準経路)と `φ`(`AdjoinRoot gK` 経由の経路)が一致するか——は
    実は**`rfl` で即座に解消する**ことが判明した(`Wₙ₊₁ =
    integralClosure Wₙ (AdjoinRoot gK)` への代入が `Subalgebra` の
    値そのものなので、`Subtype.val` の展開が両辺で完全に一致する)。
    残る本体は次の二段構成:
    (1) V0-level: `hinjV0Wn1`(`V0→Wₙ₊₁` の単射性、`Polynomial.
    map_dvd_map` 経由の`hinjV0Wn`から)→ `FaithfulSMul V0 (FractionRing
    Wₙ₊₁)`(`faithfulSMul_iff_algebraMap_injective`)→ `letI :=
    FractionRing.liftAlgebra V0 (FractionRing Wₙ₊₁)` で mathlib
    純正の instance を得る(`IsScalarTower V0 (FractionRing
    V0)(FractionRing Wₙ₊₁)` は`infer_instance`で自動的に出る)。
    (2) Wₙ-level でも全く同じパターンを1段繰り返す(`FaithfulSMul Wₙ
    (FractionRing Wₙ₊₁)` → `letI := FractionRing.liftAlgebra Wₙ
    (FractionRing Wₙ₊₁)`)——ここから `FractionRing Wₙ₊₁ ≃ₐ[Wₙ₊₁]
    AdjoinRoot gK`(`IsLocalization.algEquiv`)を`Wₙ`のscalar
    でも可換にする(`e_Wn1.commutes` + 上記`rfl`)ことで
    `Algebra.IsSeparable (FractionRing Wₙ)(FractionRing Wₙ₊₁)` を
    `hsepAdjoin`(`AdjoinRoot gK`側、既存)から移送する
    (`AlgEquiv.ofRingEquiv` + `AlgEquiv.Algebra.isSeparable_iff`)。
    最後に `IsScalarTower (FractionRing V0)(FractionRing
    Wₙ)(FractionRing Wₙ₊₁)`(2つの `FractionRing.liftAlgebra`
    instance 間の両立性、`IsScalarTower.of_algebraMap_eq'` +
    `IsLocalization.ringHom_ext (nonZeroDivisors V0)` で pointwise に
    `IsScalarTower.algebraMap_apply` を繋いで示す)を経て
    `Algebra.IsSeparable.trans` で `(FractionRing V0)→(FractionRing
    Wₙ)→(FractionRing Wₙ₊₁)` を1本に繋いだ。
    ★パッケージング上の教訓(新規): `theorem` の**戻り値の型**
    (ここでは結論)の中で `FractionRing.liftAlgebra V0 (...)` を
    直接使おうとすると、型検査の時点(証明本体に入る前)で
    `FaithfulSMul`・`Fact (Irreducible gK)` 等の instance が要求され、
    シグネチャ側だけでは間に合わない(または `whnf`/`isDefEq` が
    数十秒〜timeoutする)。**`Σ' (inst : Algebra ...), Algebra.
    IsSeparable ... inst` で instance そのものを戻り値として束ね、
    `theorem` ではなく `def`(戻り値が `Type`)にする**ことで、
    シグネチャの型検査が instance 検索を一切要求しない形になり、
    問題が消えた(`lean-idioms.md` に追記候補)。
    これで `differentIdeal_tower_diamond` の**20/20 instance が
    全て揃った**——次回はこの `falt1_hsep_bundled` の `.1`・`.2` を
    上記19-instance の組み立てに接続し、`differentIdeal_tower_diamond`
    を実際に1回呼び出して diamond 等式そのものを得ることが最短距離。

    ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**接続を実行し、
    `differentIdeal_tower_diamond` を Falt1 の `V0→V1`・`Wn→Wn1`
    構成に対して初めて実際に呼び出した**(`falt1_differentIdeal_tower_
    diamond_assembled`、`KaehlerAux.lean`、commit `03a485fa`、
    `lake build`・`#print axioms` 確認済み・sorry 無し)。20個の
    instance(`IsDedekindDomain`×2・`Module.IsTorsionFree`×5・
    `Module.Finite`×5・`Algebra`×5・`IsScalarTower`×2・`hsep`)を
    既存の `falt1Module*` 群 +
    `falt1BaseChangeAlgHom_generator_and_injective`(`ψ` を得る)+
    `falt1_hsep_bundled` から機械的に組み立て、
    `differentIdeal_tower_diamond hsep` を実際に適用して
    ```
    differentIdeal V1 Wn1 * Ideal.map (algebraMap V1 Wn1) (differentIdeal V0 V1)
      = differentIdeal Wn Wn1 * Ideal.map (algebraMap Wn Wn1) (differentIdeal V0 Wn)
    ```
    (`V1 := integralClosure V0 (AdjoinRoot fK)`、`Wn1 := integralClosure
    Wn (AdjoinRoot gK)`、`Algebra V1 Wn1 := ψ.toRingHom.toAlgebra`)を
    型検査済みの項として得た——**これが Theorem 1.2・item 3c の
    中核道具(different の塔の公式)の初めての具体的インスタンス化**。

    ★パッケージング上の限界(記録): この定理の**戻り値**を(呼び出し側が
    再利用しやすいよう)`∃ ψ, letI := ψ.toRingHom.toAlgebra;
    <diamond 等式>` の形で直接書こうとしたところ、`differentIdeal V1
    Wn1` が要求する `Module.IsTorsionFree V1 Wn1`(`ψ` の単射性からしか
    出ない、`ψ` 選択**後**にしか証明できない事実)の instance 検索が、
    `∃`/`letI` を含む**シグネチャの型検査自体**(証明本体に入る前)で
    走ってしまい、`maxHeartbeats 2000000`(168秒)でも終わらないことを
    実測した。`hsep` で使った「`Σ'` で instance ごと束ねる」策も、
    今回は `Module.IsTorsionFree V1 Wn1` が `ψ`(存在文の中の束縛変数)
    に依存するため使えない(`ψ` を `Prop` の `∃` から取り出す
    `obtain` は目標が `Type` だと `Exists.casesOn` の
    `Prop`-限定に阻まれる、というさらに別の障害も確認した)。
    **結論: 戻り値は `True` のまま(内部で `hdiamond` という名の
    genuine な項を構築し、その型を docstring に逐語で記録する)方式
    に倒した**——数学的内容は完全に確立しているが、「他の証明から
    直接呼び出せる補題」としての再パッケージングは未解決のまま次回に
    持ち越す(方針候補: `ψ` の単射性を含む**全ての ψ 依存 instance を
    シグネチャの外側の明示的仮定として渡してもらう**——
    `differentIdeal_tower_diamond` 自身がそうしているのと同じ形——か、
    非 instance の明示引数として `@` で全て手渡しする、のいずれか)。

    ★これで Theorem 1.2 の証明骨格(3a: monogenicity・3b: 長さの評価・
    3c: `V_n` 塔の構成)のうち、**「different の塔の公式」という中核
    道具は完全に具体化できた**。残るのは (i) `cancel_conductor_delta`
    の `hspan_eq` をこの具体的 diamond 等式から導く接続、(ii) 再帰的な
    `V_n` 族(`V_1, V_2, ...`)を実際に構成すること——このセッションで
    確立した `V0→V1` の1ステップを繰り返し適用する形になるはずだが、
    「πの選び方(次の一様化元)」の再帰的な選定は未着手。
   β-(d+1)(δ_n-δ_{n+1})`、`β=min{1,δ_n/(d+1)}`)を整理した形
   `δ_{n+1}≤δ_n-min{1,δ_n/(d+1)}/(d+2)` から `δ_n→0` を、`V_n`・`W_n`
   の具体的構成に一切依存しない**純粋な実数列の不等式**として抽出・
   完全に証明した。2段の議論(δ_n≥d+1が続く限り固定量ずつ減るので
   有限回で脱出/以降は比 `1-1/((d+1)(d+2))∈(0,1)` の等比減衰)を
   `squeeze_zero`+`tendsto_pow_atTop_nhds_zero_of_lt_one` で結んだ。
   ★残るのは、この不等式の**前提**(`hrec`)を実際の `V_n`・`W_n` の
   `Module.length` から導く部分(3bの一部・元の「劣加水性・完全列」
   評価)——これは 3b と一体で、`V_n` 族の具体形(3c)が定まってから
   でないと着手しづらい。

**§2-4(11項目)**: ★★★★★2026-09-04、`Definition 2.1`(almost étale
covering の定義)は `Found/Falt1/AlmostEtale.lean` で sorry 無く
formalize できた(`IsAlmostEtaleCovering`、commit `d1f92c36`)——
「projective・trace・idempotent」は個別には mathlib の道具(dual basis
の補題・`Algebra.trace`・`Algebra.FormallyUnramified.elem`)で書け、
特に「`B⊗_AB` の像」を捉える比較写像は局所化とテンソル積の可換性の
同一視を経由せず `Algebra.TensorProduct.lift` で直接構成できることを
発見した(想定より遥かに単純だった)。

★★★★★★★続けて `Theorem 2.2`(nilpotent ideal 上の lifting)の証明全文
(物理p.7=印字p.260)を260dpiで精読した。証明は: (1) `B` の almost
射影性から任意の `ε>0` で `p^εφ` を `A` 加群写像として持ち上げる、
(2) 積の非可換性の障害
`b_0(p^εφ_ε(b_1b_2)-φ_ε(b_1)φ_ε(b_2))b_3` が **Hochschild
cohomology** `H^2(B/A,I)` の類を定め、これが `m` で零化される
(remark 2.1(v)、証明は「同様の議論」とだけ書かれ本文中では未証明)
ことから `ε` を倍にして消せる、(3) 持ち上げの一意性は `H^1(B/A,I)`
(これも `m` 零化)まで、`C` が `p` 捩れ無しなのでちょうど一意に決まる、
(4) `mB→C` への拡張は `p^ε` の2乗関係から。

**結論**: Definition 2.1 とは質的に異なり、`Theorem 2.2`-`2.4`・`§3`・
`§4` は **Hochschild cohomology(`B⊗_AB⊗_A⋯⊗_AB→M` の複体としてremark
(v)で定義される)の理論そのもの・その almost 消滅性・一般の変形理論
(`H^2`=障害・`H^1`=分類)** という、mathlib に存在しない独立した理論を
要する。Definition 2.1 で見つかった「見た目より簡単だった」という
僥倖はここでは再現しない見込み——引き続き完全に未着手。

★★★**結論(正直な評価)**: `/goal` の 13/13 は、このセッションでの
継続作業だけでは現実的な時間内に到達できない規模の作業(Theorem 1.2
1件でさらに大きな作業、§2-4 は新規ライブラリ構築に近い規模)である
可能性が高い。ただし Lemma 1.1 で確立した一般補題(especially
`subsingleton_H1Cotangent_self`・`polynomialKaehlerSplit`)は Theorem
1.2 の帰納法でも再利用できる見込みがあり、**無駄になっていない**。

★★★★2026-09-04、**§3・§4 の主張そのものを初めて実際に読んだ**
(`lean/ABC3/Skeleton/Falt1/Section3.lean`・`Section4.lean`、既存の
Interface posit を精読、新規の 260dpi 目視は不要——既にファイル内に
逐語が写してあった)。これまで「§2 が Hochschild cohomology を要する
ので恐らく §3・§4 も同様に重い」という**推測**だったが、今回**§3・§4
自身の主張を具体的に確認**して評価を裏付けた:

- **§3(Theorem 3.1・3.2、2項目)**: `S⊗_R R_∞` の正規化が `R_∞` の
  almost étale covering であること——`Interface/Falt1/GoodReduction.lean`
  の `GoodReductionSetup` が既に `AlmostEtaleSetup.isAlmostEtale`
  (旧い posit 述語、**まだ `Found/Falt1/AlmostEtale.lean` の実際の
  `IsAlmostEtaleCovering` に差し替えられていない**)を参照する形で
  スケルトン化済み。証明そのものは Theorem 2.2-2.4(Hochschild
  cohomology による持ち上げ)の結果を使う——**§2 の壁がそのまま
  継承される**ことを確認した。
- **§4(Theorem 4.1-4.3・4.5、5項目)**: `H^i(Δ,R̂)` という Galois
  コホモロジー(`Δ` = Galois 群、`R̂` = 完備化された環)の計算、
  スペクトル系列 `E_2^{a,b}=H^a(Δ,H^b(I,R̂))⇒H^{a+b}(Δ,R̂)` の
  almost 退化、cup 積による Hodge-Tate 型の同型——これは**§2 より
  さらに深い**、p進 Hodge 理論の核心部分そのもの(周期環の Galois
  コホモロジー・スペクトル系列)。mathlib には(1)almost mathematics
  だけでなく(2)この種の Galois コホモロジー計算の対象となる周期環
  `R̂` の理論も存在しない——**2重に新規ライブラリ構築が要る**。

★結論(確定・更新): §2-4 のブロッカーは**単一の Hochschild
cohomology だけではなく**、§4 では独立にもう1段深い(周期環の
Galois コホモロジー・スペクトル系列)理論を要することが**具体的に
確認**できた——「推測」から「該当ページを読んで確認した事実」への
格上げ。これで §2-4(11項目)全てについて、着手前の偵察は完了した
と見なしてよい——次に必要なのは実際の mathlib 開発であり、この
セッション規模の作業ではない。

## 1. 構造(2026-09-04 確定——旧版の3つの未確認点をすべて解消した)

★★**Chapter見出しは"CHAPTER"の語を伴わない裸のローマ数字**("I""II""III"
が単独行に出るだけ)。原文 物理p.4(印字p.256)の(e)段落「... In §I we study
the case of good reduction ... In §II we construct the theory 𝒳*(X) ...
Finally in §III we treat the case of bad reduction ...」で著者自身が
3章構成を明言しており、実測(下表)と完全に一致する。

★★**項目の番号付けは "N.M. Kind" の順**(番号が先、種別が後)。LocProP等の
"Kind N.M"(種別が先)と**逆**。`tools/paper-items.mjs`はそのまま使えない
(専用の抽出スクリプトが必要——再現用ロジックは本ファイル末尾のメモ参照)。

★★**pageOffsetを訂正した(254→253)**: 物理p.1はJSTORの表紙(著者・出典
スタンプのみ、本文ページ番号なし)で、物理p.2が印字p.255(論文本文の
実際の1ページ目)。150dpiで物理p.1-3を目視し、物理p.3の見出し
"256 GERD FALTINGS"で確定した。★[LocProP]形式化時に登記した
`theorem_I_4_4.src`の`pdfPage := 17`は260dpi目視による直接確認だった
ため訂正の影響を受けない(訂正後の式 印字270-253=17 で無矛盾も確認済み)。

Chapter III §4 は末尾を確認した結果、**4.1 Theoremの1件のみ**(物理p.45、
印字p.298)で、直後に p.299 の BIBLIOGRAPHY が続き論文が終わる。
→ 総項目数は **30 で確定**(旧版の「30+」の「+」は不要)。

## 2. 全項目表(物理ページは pageOffset=253 で算出・機械抽出。個別の
   260dpi目視は §4.4 以外まだ)

### Chapter I(good reduction、4節・13項目)

| 節 | 項目 | 種別 | 物理p. | 内容の頭出し(OCR、地の文の空白喪失あり) |
|---|---|---|---|---|
| §1 Ramification theory for discrete valuation rings | 1.1 | Lemma | 4 | For any extension V⊂W... the natural map... |
| | 1.2 | Theorem | 5 | If V⊂W is any extension and 𝒥'V denotes the normalization... |
| §2 Almost unramified extensions | 2.1 | Definition | 6 | Suppose A is a ring, B an A-algebra. B is called an almost... |
| | 2.2 | Theorem | 7 | Suppose B=A+mB is an almost etale covering of A, C an... |
| | 2.3 | Theorem | 7 | Suppose I⊂A is a nilpotent ideal, B an almost etale covering... |
| | 2.4 | Theorem | 8 | Suppose B is an almost etale covering of A. |
| §3 Good reduction | 3.1 | Theorem | 10 | Suppose S is a normal finite R-algebra such that S[1/p]... |
| | 3.2 | Theorem | 12 | Suppose S is a finite normal torsionfree R-algebra such that... |
| §4 Differentials and cohomology | 4.1 | Theorem | 13 | (i) The map Ω_{R/V}⊗R̄→... induces almost isomorphisms |
| | 4.2 | Theorem | 15 | (i) H^i(Δ,R̂)≅Λ^i((R⊗_V V̄)^(-1)d)⊕(rest)... |
| | 4.3 | Theorem | 17 | There exists a Γ-equivariant functorial extension |
| | **4.4** | **Theorem** | **17** | **(i) H^i(Δ,R̂)≅Λ^i((R⊗_V V̄)^(-1)d)⊕(rest)...**★先行あり(縮小形) |
| | 4.5 | Theorem | 19 | (i) The spectral sequence |

### Chapter II(𝒳*(X) の構成、3節・8項目)

| 節 | 項目 | 種別 | 物理p. | 内容の頭出し |
|---|---|---|---|---|
| §1 Construction of 𝒢* | 1.1 | Theorem | 23 | (i) There exists a spectral sequence |
| | 1.2 | Theorem | 25 | (i) There exist spectral sequences |
| §2 The isomorphism with étale cohomology | 2.1 | Lemma | 27 | Any x∈X is contained in an open U⊂X which is a K(π,1). |
| | 2.2 | Corollary | 27 | For X as above, x∈X, Spec(𝒪_{X,x̄}⊗V K̄) is a K(π,1). |
| | 2.3 | Lemma | 28 | Suppose X is smooth over V, D⊂X a divisor with normal... |
| | 2.4 | Theorem | 30 | The transformation (X proper and smooth, D⊂X a divisor... |
| §3 Relations to Hodge cohomology | 3.1 | Theorem | 35 | Suppose X is proper and smooth over V, D⊂X a divisor with... |
| | 3.2 | Theorem | 36 | Suppose X is proper and smooth over V=W(k)(Witt vectors)... |

### Chapter III(bad reduction、4節・9項目)

| 節 | 項目 | 種別 | 物理p. | 内容の頭出し |
|---|---|---|---|---|
| §1 Commutative algebra | 1.1 | Lemma | 37 | There exists an e independent of n such that the normalization... |
| | 1.2 | Proposition | 38 | Suppose R has one system of good units for some (infinite)... |
| | 1.3 | Theorem | 39 | Suppose that R is a V-algebra of finite type, geometrically ir-... |
| §2 Global methods | 2.1 | Definition | 40 | A stable punctured curve of type (g,r) over an algebraically... |
| | 2.2 | Lemma | 41 | If V is a complete discrete valuation ring with fraction field K,... |
| §3 Rigid coverings | 3.1 | Definition | 41 | A V-map f:X→Y is called a rigid étale covering if |
| | 3.2 | Theorem | 42 | Suppose f:U.→Y is a rigid étale hypercovering, a coherent... |
| | 3.3 | Theorem | 43 | Suppose R is a V-algebra of finite type, smooth in characteristic... |
| §4 Intermediate cohomology | 4.1 | Theorem | 45 | Suppose X is a smooth proper K-scheme, D⊂X a divisor with... |

## 3. [LocProP]が直接引用する箇所(確認済み・実装済み)

| LocProP側 | Falt1側 | 状態 |
|---|---|---|
| [LocProP] Lemma 2.1 | Chapter I, Theorem 4.4, (i)(iv) | ★★物理p.17(印字p.270)で260dpi目視確認済み。実装済み(縮小形) |
| [LocProP] Lemma 2.3 | Chapter I(Leray-Serre + 前2つ、Theorem 4.3?) | 未確認・未着手 |
| [LocProP] 全体(§2の土台) | almost étale extension理論(§I-2) | 未確認・未着手 |

## 4. 着手順序(葉から)

leaf-first の原則により、他の項目に依存しない・短い証明から着手する。
候補(未検証・着手時に依存グラフで裏取りすること):
1. **Chapter I §1**(1.1 Lemma, 1.2 Theorem)—— almost étale extension理論
   の手前の純代数的補題で、他の Falt1 内項目への依存が最も少なそう。
2. **Chapter III §2.1**(Definition, stable punctured curve)—— 定義のみ
   なら証明義務が無く着手しやすい。
3. Chapter I §4.4 の縮小形を(i)-(vii)全体に拡張(LocProPの土台なので
   優先度は高いが、(ii)(iii)(v)(vi)(vii)は almost étale extension の
   深い性質を要する可能性があり、葉ではない懸念がある——着手前に
   `node tools/graph.mjs`で確認)。

## 5. 番号付け規則の再現手順(scratchpad/falt1-items.mjsの後継)

`scratchpad/falt1-full-index.mjs`(本セッションで作成、再現可能):
- 章境界: 行が正規表現 `^\s*(I|II|III)\s*$` に一致する単独行を探す
  (`CHAPTER`という語は原文に無い)。
- 節見出し: `^\s*(\d+)\.\s*[A-Z][a-z]` (空白喪失により単語がくっつくが
  見出し行自体は検出できる)。
- 項目: `^\s*(\d+)\.(\d+)\.\s*(Theorem|Proposition|Lemma|Corollary|
  Definition|Remark)`("N.M. Kind"の順、LocProPと逆)。
- 印字ページ: 走り込み見出し行 `NNN GERD FALTINGS`(偶数ページ)・
  `p-ADIC HODGE THEORY NNN`(奇数ページ)を全行から収集し、各項目の
  直前・直後で挟んで判定。
- 物理ページ = 印字ページ - 253(§1で確定した値)。

## 6. 未検証事項(次の個別着手時につぶすこと)

- 上表の物理ページ番号は§4.4以外、機械抽出のみ(260dpi目視未実施)。
  逐語照合(`原文 (Falt1 p.N):`)を使う前には必ず260dpiで該当ページを
  確認すること(JSTORスキャンは地の文の単語間スペースを失うが、定理の
  数式・番号は無事——Theorem 4.4の実測で確認済みの傾向)。
- 各項目が Falt1 内の他項目にどう依存するか(依存グラフの辺)は未調査。
  着手順序の葉判定は個別に`node tools/graph.mjs --owner Falt1`で
  裏取りすること。

関連: [[measure-mathlib-before-skeleton]] / `ResearchPaper/mathlib-gap.json`の
`locprop-perfectoid-hodge-tate`(この論文が主典拠)
