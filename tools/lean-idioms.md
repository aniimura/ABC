# 配管の対策集 —— この codebase で繰り返し出る失敗形と、効いた直し方

★これは**数学の難しさではなく、エラボレータとの戦い方**の記録である。
同じ穴に何度も落ちたので、1 箇所にまとめる。

★★**運用**: 新しい失敗形に当たったら、直すと同時にここへ 1 行足すこと。
「前にも見た」と思ったらまずここを引く。

---

## 1. `instances 透明度で型が合わない` —— 最頻出

**症状**: `rw` / `simp` / インスタンス探索が
「`Did not find an occurrence of the pattern ...`」と言って落ちる。
目標にはその式が**そのまま見えている**。エラー末尾に

```
Note: The target expression is not type-correct under the `instances` transparency level
```

が付く。原因は `Under.right W` / `Z.unop.left.obj` / `WideSubcategory` /
`PfRootObj` の射影が、`instances` 透明度では展開されず型が合わなくなること。

**効いた直し方(上から順に試す)**

| 手 | 例 |
|---|---|
| ★★**まずこれを試す**: `set_option backward.isDefEq.respectTransparency false in` を宣言の直前に置く | mathlib の `AffineTransitionLimit.lean` 自身が同種の場面(`Subtype`包みの`def`をPreorder-as-categoryの対象に使う等)で多用している。`rw`/`simp`の`Category.comp_id`等が「unused」と言われるのに式が消えない、という症状にも効く(`CorrHyp/FieldLimit.lean::cocone_ι_congr`、2026-09-04)。**副作用**: 宣言全体でdefeq判定が緩むので、他の箇所で意図しない項が通ってしまわないか要確認——`lake build`が通ることと`#print axioms`がsorryAxを含まないことは必ず確認する。 |
| `Matrix.SpecialLinearGroup n R := {A // A.det=1}` の匿名コンストラクタ`⟨N,hNdet⟩`が期待した`SpecialLinearGroup`型ではなく素の`Subtype`として elaborate される(`Application type mismatch`) | ★★**解決した(2026-09-04、`CorrHyp` Hecke型共役の計算で発見)**: `M : Matrix.SpecialLinearGroup n R`が手元にあるとき、`(M : Matrix.SpecialLinearGroup n R) i j`のように**`M`自身をその場で何度も添字アクセスしない**——最初に`set A : Matrix n n R := (M : Matrix n n R) with hA`で台の`Matrix`を1回だけ取り出し、以後は`A`だけを使って計算する(`Gamma_mem`等から得た性質も`rw [← hA] at ...`で`A`の言葉に変換してから使う)。証明の最後、`SpecialLinearGroup`の等式を示す段でも`apply Matrix.SpecialLinearGroup.ext; intro i; intro j; fin_cases i <;> fin_cases j <;> simp_all [...]`のように**`Matrix`レベルの等式に早く落として**そこで戦う。`Matrix.SpecialLinearGroup.coe_pow`等の coe補題を`rw`で当てるときも、対象の式が`↑(g^k)`(先にpow、後でcoe)なのか`(↑g)^k`(先にcoe)なのか型注釈の書き方でLeanの中間elaborationが揺れる——`#check`で実際の形を確認してから、必要なら`((g^k : SpecialLinearGroup n R) : Matrix n n R) = ...`のように**多重に型注釈して意図した括り方を強制する**。 |
| ★`rw` をやめて **`Eq.trans` / `congrArg` の項**で繋ぐ | `(Category.assoc _ _ _).trans (h.trans (Category.assoc _ _ _).symm)` |
| ★`congrArg` には**関数の型を明示**する | `congrArg (fun t : X ⟶ Y => t ≫ f) h` |
| ★`show` で**きれいな型**に言い換えてから触る | `show HomBirat.mk (biratPfIdx …) (rootMap … ≫ …) = _` |
| 射影を**名前つき `def`** に逃がす | `rtObjPf` / `pfDown`(既存の実例) |
| 添字対象を `idxBiratMk P G a hc hs` の形で**受け取る** | 型が素の `C` の対象になる |

## 2. `cancel_epi` / `cancel_mono` がインスタンスを見つけない

**症状**: `have hep : Epi f := …` を置いたのに
`(cancel_epi _).mp h` が `failed to synthesize Epi f` で落ちる
(`haveI` にしても同じ)。

**直し方**: **構造体の射影を直接使う**。

```lean
have hep : Epi f := P.totEpiC _ _ f
exact hep.left_cancellation _ _ h        -- cancel_epi の代わり
exact hmono.right_cancellation _ _ h     -- cancel_mono の代わり
```

★`Epi` / `Mono` は `left_cancellation` / `right_cancellation` を
フィールドに持つので、インスタンス探索を経由せずに済む。

## 3. `ℕ≥1` が**依存位置**(対象の中)に現れる

**症状**: `⟨A, n⟩ : PfRootObj P F` のように次数が対象の一部になっていると、
「次数は命題として等しいが項としては違う」2 つの構成が**同じ項にならない**。
`rw` は `motive is not type correct` で落ちる。

**直し方**: ★★**次数を仮引数 `K` に出し、等式を仮引数で受け取る**。

```lean
def fooK (W : Idx …) (K : ℕ+) (hK : degFr W = K) : … := …   -- 一般形
theorem fooK_eq … : fooK W K hK … = foo W … := by subst hK; rfl
```

`K` は**変数**なので `subst` が効き、証明部分は `Prop` なので
**証明が違っても同じ項**になる。
実例: `biratPfIsoA'` / `biratPfMk'` / `biratPfMk'_eq`(`Prop55BiratPf.lean`)。

★**規約**: 依存位置に `ℕ≥1` を持つ定義を書いたら、
**その場で primed 版(次数を仮引数に出した版)も並置する**。
後から足すと呼び出し側を全部書き直すことになる。

## 4. `subst` を使うために、等式の片側を**変数**にする

**症状**: `h : f x = c.field` の `c.field` は変数ではないので `subst` できない。

**直し方**: 補題の仮引数に `δ` を取り、`hδ : f x = δ` の形で受ける。
呼び出し側で `δ := c.field` を渡せばよい。
実例: `biratPfHom_surj_mk`(`δ` と `φ` を仮引数に出してから `subst`)。

★構造体の **eta** が効くので `idxBiratMk P G T.unop.hom.hom _ _ = T` は **`rfl`** である
(`Over` / `WideSubcategory` / `Discrete PUnit` / `Opposite` すべて eta を持つ)。
これで「一般の添字対象」を「構成子の形」に**無料で**言い換えられる。

## 5. `rw` の末尾 `rfl` が閉じない

**症状**: 目標が見た目 `X * Y = X * Y` なのに `unsolved goals`。

**直し方**: `rw` の自動 `rfl` は `reducible` 透明度なので、
`exact rfl` / `exact mul_comm _ _` / `exact congrArg (· * y) h` と**項で閉じる**。

## 6. `𝟙` がどの圏の恒等射か決まらない

**症状**: `PfCat P F` は定義上 `C` そのものなので、
`𝟙 (X : PfCat P F)` と書いてもインスタンスは `C` の方が選ばれる。
型注釈を対象に付けても直らない。

**直し方**: **恒等射を像として書く** —— `toHomPf (𝟙 X)` / `toHomBirat (𝟙 X)`。
関手の `map_id` が `rfl` なので、これは本当に恒等射である。
実例: `biratPfHom_id`。

## 7. 巨大な構造体を型の中に直接書かない

**症状**: `Frobenioid (pfRootPre P F)` を型の中に書くと
`whnf` が 200000 heartbeats を超える(**実測 73 秒で timeout**)。

**直し方**: **仮引数で受ける**(`(Gpf : Frobenioid (pfRootPre P F))`)。0.2 秒になる。

## 8. `Under.isoMk` / 無名構成子が未簡約の型に当たらない

**直し方**: 同型を**型注釈つきの `let`** で別に組んでから使う。

```lean
let e : W' ≅ d.right := { hom := …, inv := …, … }
```

## 9. 宇宙変数の使い回し

**症状**: 始域と終域で別の圏を扱う補題で `u2 v2` を使い回すと当たらない。

**直し方**: `universe uu1 vv1 uu2 vv2` と**新しく宣言する**。
実例: `iso_unique_of_rigid`(`Cor54Rigid.lean`)。

## 10. ファイル書き込みの事故

* Lean のファイル内容は **`Write` ツールで書く**。
  python の heredoc に `𝒞` のようなサロゲート対のエスケープを入れると
  Windows の既定エンコーディングで `UnicodeEncodeError` になり、
  **ファイルが途中まで書かれて壊れる**(実際に 772 行失った)。
* `git commit -m` にバッククォートを入れない。`git commit -F -` ＋ heredoc を使う。
* `lake build` は `lean/` から、`node tools/check.mjs` は**リポジトリ直下**から。

## 11. doc コメントと宣言の間に `set_option ... in` を置けない

**症状**: `unexpected token 'set_option'; expected 'lemma'`。
MCP(`lean_check`)では通るのに `lake build` で落ちるのでたちが悪い
——MCP へ投げる断片には doc コメントを付けないことが多いため。

```lean
/-- ★説明 -/
set_option maxHeartbeats 1000000 in
theorem foo : … := …          -- ✗ パースエラー
```

**直し方**: `set_option ... in` を **doc コメントより前**に書く。

```lean
set_option maxHeartbeats 1000000 in
/-- ★説明 -/
theorem foo : … := …          -- ✓
```

`variable (P F) in` / `omit … in` も同じ。★**属性 `@[simp]` は逆で、doc の後**である。

## 12. 「包み `def` の射影」が `rw` の照合を止める

**症状**: 補題は `X.obj` について述べているのに、目標には
`(scaleRootObj k X).obj` と出ていて `rw` が当たらない。定義上は同じ
(`scaleRootObj k X := ⟨X.obj, k * X.root⟩`)だが**構文が違う**。
1 と同じ `instances` 透明度の注記が付くことが多い。

**直し方**: 触る直前に**包みだけ `unfold`** する。

```lean
rw [hR, hL]
unfold scaleRootHom scaleRootObj   -- ← これで (scaleRootObj k X).obj が X.obj になる
rw [hf, hg]
```

★★`simp only [scaleRootObj]` でも同じだが、`unfold` の方が
「包みを 1 枚剥いだだけ」であることが読み手に分かる。
★依存位置の `ℕ≥1`(3 番)と合わせて、**包みの中の根は必ず項として書き下せる形にしておく**。

## 13. `IsIso` のインスタンスが `haveI` で登録されない

**症状**: 仮定に `ha : IsIso (P.Base a)` があり `haveI := ha` もしたのに、
`IsIso.hom_inv_id` などが `failed to synthesize instance IsIso (P.Base a)` と言う。
★名前をつけて `haveI hA : IsIso (P.Base a) := ha` としても直らないことがある
(`rw` で式を書き換えた**後**の目標では、射の項が `instances` 透明度で
別物に見えているため)。

**直し方**: **インスタンスを明示的に渡す**。

```lean
exact @IsIso.hom_inv_id _ _ _ _ (P.Base a) ha       -- ✓ 確実
-- exact IsIso.hom_inv_id _                          -- ✗ 合成に失敗しうる
```

★`inv` を含む等式を述べるときも `@inv _ _ _ _ f ha` と書いておくと、
呼び出し側でインスタンスがずれない。`IsIso` は `Prop` なので
**証明無関係により、どの証明を渡しても同じ項**になる。

---

## 検査器のキャッシュ(2026-08-21)

`node tools/check.mjs` の時間はほぼ全部が `pdftotext` の呼び出しだった。
`.cache/pdf-pages.json` に**ページ本文を跨いでキャッシュ**するようにした
(**55 秒 → 7 秒**)。

★鍵は `PDF のパス # ページ # mtime # size` に加えて
**`check.mjs` 自身のハッシュ**を含む —— `normalize` / `squash` / `PDF_MODES` を
触ったら必ず外れる。これを忘れると
「正規化を変えたのに古いテキストで通る」という**器具の穴**になる。

---

## 14. 型の同義語(`PfCat P F := C`)を跨ぐと `rw` が当たらない(2026-08-24)

**失敗形**: `PfCat P F` は `C` の型の同義語で、`pfDiv` などは
暗黙引数を `{A B : PfCat P F}` で取る。一方 `toHomPf ψ` は `ψ : A ⟶ B`
(`A B : C`)から作るので、**同じ項でも暗黙引数の書かれ方が違う**。

```
-- 目標側:   pfDiv (A := pfDown P F ((toPfCat P F).obj X)) …
-- 補題側:   pfDiv (A := rtObj P F (pfObjDown P F X) 1) …
```

`rw` は `instances` 透明度でしか合わせないので、
**defeq なのに「パターンが見つからない」**になる。

★`(toPfCat P F).obj A = A` の simp 補題を足しても、今度は逆向きに
「`C` の項が `PfCat P F` を期待されている」でずれる。

**対処(測定済み)**:
* 同義語の**両向き**に名前付きの `abbrev` を置く(`pfObjDown` / `pfUp`)。
* それでも駄目なら `rw` をやめて `Eq.trans` / `congrArg` の**項の側**で組む
  (`exact` は既定透明度で合わせるので通る)。
* 圏の合成 `≫` を `compPf` に開く橋 `pfCat_comp_eq : f ≫ g = compPf P F f g := rfl`
  を用意しておくと `pfDiv_comp` 等が当たるようになる(ただし上のずれは別途残る)。

★**2026-08-24 の未着手**: `pf_isOfIsotropicType`(`𝒞^pf` の根 1 の部分が
isotropic 型)はこの穴で止まっている。段取りは
`Found/FrdI/Prop53PfCatRoot.lean` の ★2 に書いてある。

### 14 の続報(2026-08-24、2 回目の試作)

**主因は「自前の同義語ほどき関数を書いたこと」**だった。
在庫に `pfDown (A : PfCat P F) : C := A` があるのに `pfObjDown` を自分で足すと、
`pfDiv` などの暗黙引数(在庫の `pfDown` で書かれている)と噛み合わない。

★**在庫の同義語ほどき関数に揃える**。揃えた上で

* `pfCat_comp_eq : f ≫ g = compPf P F f g := rfl`
* `toPfCat_map_eq : (toPfCat P F).map ψ = toHomPf ψ := rfl`

の 2 本を橋にすると、`pfDiv_comp` / `pfDeg_comp` / `rootDiv` / `rootBase` /
`compRoot_root_one` はすべて当たるようになった。

★★**目標の側で `rw` しない** —— `≫` を開いた目標は
`instances` 透明度で型が付かなくなる(エラーに
「The target expression is not type-correct under the `instances` transparency level」
と出る)。**自分で述べた `have` の側で `rw` して、最後に `exact` で defeq に頼る。**

★残るのは `IsIso` のインスタンス合成(既存の仮定が暗黙引数の書かれ方の違いで
拾われない)で、これは `haveI` を**目標に現れる形そのまま**で置き直すしかない。

### 14 の続報(2026-08-25、3 回目)—— `𝒟` 側の合成が `Functor.comp` をまたぐとき

`Cor54SeamCls.lean` で、**底 `𝒟` の射の合成**でも同じ穴に落ちた。
`sq : P₁.proj ⋙ ΨB ≅ Ψ ⋙ P₂.proj` の成分 `sq.hom.app X` の終域は
`(Ψ ⋙ P₂.proj).obj X` だが、`P₂.Base f` の始域は `(P₂.toElem.obj (Ψ.obj X)).base`。
**defeq だが `instances` 透明度では型が付かない**ので、

* `Category.assoc` / `IsIso.hom_inv_id` が `rw` で当たらない
* `IsIso (P₂.Base (Ψ.map φ))` が `haveI` で置いてあっても**合成に失敗する**

**対処(測定済み、3 点)**:

1. **`inv` は `@inv _ _ _ _ f h` でインスタンスを明示**する。
   `IsIso` は `Prop` なので、どの証明を渡しても defeq。合成に頼らないのが確実。
2. **`Category.assoc` / `hom_inv_id` は `have` で項として置く**
   (`(Category.assoc _ _ _).symm` / `@IsIso.hom_inv_id _ _ _ _ f h`)。
   置いた `have` を `rw` するか、`Eq.trans` / `congrArg` で項の側から組む。
3. **`rw` の末尾で閉じない**ことがある —— 目標が見た目 `X = X` になっても
   `rw` の自動 `rfl` は reducible 透明度なので通らない。**`rfl` を 1 行足す**。

★★**シェルの罠(同日)**: `perl -0pi -e 's/.../(@inv ...)/'` の置換文字列で
`@inv` が**配列展開されて消える**。置換に `@` を含めるときは `sed` を使うか
`\@` でエスケープする。Lean ファイルの中身は Write/Edit で書くのが安全。

## 14 の続報(2026-08-25、4 回目)—— `inv` の `rw` は**抽象補題へ逃がす**のが唯一安定

`Prop44Univ.lean`(`𝒞^birat` の普遍性)で 4 つ新しい失敗形に当たった。
**3 回目までの対処(`@inv` で明示、`have` で項にする)では足りない場面がある。**

| 失敗形 | 出るメッセージ | 直し方 |
|---|---|---|
| `rw [Category.assoc]` が当たらない | `Did not find an occurrence of the pattern (?f ≫ ?g) ≫ ?h` (**目に見えて在るのに**) | ↓ |
| `rw [IsIso.hom_inv_id]` が当たらない | `Did not find an occurrence of the pattern inv ?f ≫ ?f` | ↓ |
| `rw [Ω.map_comp]` / `rw [Ω.map_id]` | `motive is not type correct`(`IsIso _a` に型が付かない) | ↓ |

★**原因**は共通で、`IsIso` の引数がメタ変数だと**インスタンス探索が走らない**こと。
`rw` は補題 `IsIso.hom_inv_id` を使うのに `[IsIso ?f]` を解決できず、
パターンそのものを作れないので「見つからない」と言う。

### ★★対処 —— **圏の中だけの抽象補題**に逃がす

`IsIso` を**インスタンス束縛**にした補題を別に立てて、`exact` で当てる。
補題の中では `?f` が本物の変数なので `rw` が普通に動く。

```lean
theorem frac_key_aux {X Y T U : E} (g : X ⟶ Y) [IsIso g] (a : X ⟶ T) (p : Y ⟶ U)
    (w : T ⟶ U) [IsIso w] (hsq : g ≫ p = a ≫ w) : inv g ≫ a = p ≫ inv w := by
  rw [IsIso.inv_comp_eq, ← Category.assoc, hsq, Category.assoc, IsIso.hom_inv_id,
    Category.comp_id]
```

**`rw [Ω.map_comp]` の motive 問題**も、等式を仮定に外出しすれば消える:

```lean
theorem frac_comp_aux ... (gz : X ⟶ Z) [IsIso gz] (aq : X ⟶ V)
    (hgz : gz = g ≫ z) (haq : aq = a ≫ q) ... := by
  subst hgz; subst haq   -- ★ここで初めて合成の形になる。inv の下を rw しない
```

### ★注意 1 —— 抽象補題の **universe を 1 本に固定する**

`D` / `C` を使わない補題を同じ section に置くと、`Category.{max v u2 v2} E` の
`max v u2 v2` が**そのまま universe 変数 3 本として一般化**され、
呼び出し側で `stuck at solving universe constraint` になる。

```lean
section FracAux
universe uv
variable {E : Type uE} [Category.{uv} E]   -- ★ 1 本にする
```

### ★注意 2 —— 呼ぶときは **インスタンス引数を `@` で明示**する

`exact frac_comp_aux g z a p w q …` は
`failed to synthesize IsIso (Ω.map …)` で落ちる。
**局所インスタンスは在る**(`have test : IsIso … := inferInstance` は通る)のに、
暗黙の**対象**がメタ変数のうちに探索が走るため。

```lean
exact @frac_comp_aux _ _ _ _ _ _ _ _ g hγ z hZ a p w hW q gz hP aq
  (Ω.map_comp _ _) (Ω.map_comp _ _) key
```

## 15. `calc` は defeq を渡らない —— `have` ＋ `Eq.trans` に置き換える

`Prop55BiratOmega.lean` で繰り返し当たった形。

```
error: invalid 'calc' step, failed to synthesize `Trans` instance
  Trans Eq Eq ?m.848
```

★**原因**: `calc` の各段は**型が構文的に一致**していないと `Trans` が付かない。
`omegaObj F F' ⟨A,k⟩` と `⟨biratUp A, k⟩` のように **defeq だが構文が違う**と落ちる。

★**逃げ方**: 段を `have s1 := …` で置き(型注釈を書かない)、
最後に `(((s1.trans s2).trans s3) … )` で繋ぐ。
`Eq.trans` は defeq で通るので、これだけで直る。

```lean
  have s1 := congrArg (compRoot Q F' _) (omegaMap_pfKappa F F' B k).symm
  have s2 := (omegaMap_comp F F' _ _).symm
  have s3 := congrArg (omegaMap F F') (rootMap_spec (F := F) hfi f k)
  exact rootMap_ext (F := F') hfiB _ _ (((s1.trans s2).trans s3).trans hR.symm)
```

## 16. `IsIso.inv_eq_of_hom_inv_id` は `@` で `f` とインスタンスを明示する

```
error: failed to synthesize instance of type class IsIso (rtExt (biratPre P G) F' (biratUp P G A) 1)
```

局所インスタンスは在るのに落ちる。★**暗黙の `f` がメタ変数のうちに探索が走る**ため
(idiom 14 の注意 2 と同じ)。

★**引数の順は `{f} [IsIso f] {g}`** —— インスタンスは `g` より**前**である。

```lean
refine (@IsIso.inv_eq_of_hom_inv_id _ _ _ _
  (rtExt (biratPre P G) F' (biratUp P G A) 1) hq1   -- ★f, インスタンスの順
  (biratRtIso F F' A 1 ≫ (toBiratCat P G).map (rtOneInv A)) ?_).symm
```

## 17. `set_option … in` は docstring の**前**に置く

`variable (F) in` と併用するときは

```lean
set_option maxHeartbeats 1600000 in
variable (F) in
/-- doc -/
theorem foo …
```

★docstring の中に紛れ込ませると本文の一部になって効かない
(`sed` で行番号を数えて挿入すると起きやすい)。

## 18. `rw` は「対象が関手の像の形」だと当たらない —— 項スタイルへ

`Prop55BiratOmega.lean` の組み立てで**丸一日ぶんの試行を溶かした形**。測定結果:

補題の側の対象が `omegaObj F F' X`(関手の対象写像を当てた形)で、
目標の側では既に `⟨biratUp Z₀, k⟩` に簡約されているとき、`rw` は

```
omegaObj F F' ?X =?= ⟨biratUp Z₀, k⟩
```

を解こうとする。`omegaObj` は構造体リテラルを返す `def` なので、
これは `?X.obj =?= Z₀` という**メタ変数の射影**に化けて解けない。
★エラーは「`Did not find an occurrence of the pattern`」と出るが、**原因はこれ**。

★★**対処**: 組み立ては `Eq.trans` / `congrArg` の**項スタイル**で書く。
`≫` 版の補題を用意しても同じ理由で `rw` は当たらない(実測)。

### ★付随して測ったこと

* `HomRoot P F X Y` は `≫` の**左**に置けば `Quiver.Hom` に解けるが、
  等式の**右辺**に単独で置くと解けない(メタ変数の postpone)。
  `≫` 版の補題を書くときは型注釈
  `(… : omegaObj F F' X ⟶ omegaObj F F' Zz)` が要る。
* `Mono ((Ω).mapIso i).hom` の**インスタンス探索は heartbeat を使い切る**。
  `Iso.cancel_iso_hom_right` を使うこと。
* 文の中で `show T from e` を使うと `have this := e; this` に展開され、
  以後の `Category.assoc` などの単一化を**全部止める**。
  文では名前つき引数 `(X := …) (Y := …)` を使い、`show` は**戦術の側**に置く。

## `letI := algOfHom f` で `Algebra L L` を入れても `Algebra (𝓞 L) (𝓞 L)` は f を見ない

**失敗形**

```lean
letI := algOfHom f   -- Algebra L.toIF L.toIF
have h : ((algebraMap (𝓞 L.toIF) (𝓞 L.toIF)) x : L.toIF) = (FinSub.hom f) (x : L.toIF) := rfl
-- Type mismatch: rfl has type ?m = ?m
```

★★原因: `Algebra (𝓞 A) (𝓞 A)` には **`Algebra.id` が既にインスタンスとしてある**ので、
`letI` で入れた `Algebra L L` から**派生する**インスタンスより先に選ばれる。
その結果 `algebraMap (𝓞 L) (𝓞 L)` は**恒等写像**に解決され、`f` を一切見ない。

★★★**ただし既存の定義（`resHOS` など）は無事**である ——
`resHOS {L M : Type} [Algebra L M]` の本体は `L`・`M` が**相異なる型変数**の状態で
エラボレートされており、そこでは `Algebra.id` は候補にならないので
**派生インスタンスが焼き付いている**。あとから `L := M := L.toIF` を代入しても
型クラス探索は**やり直されない**ので、`resHOS` は正しく `f` に依存する。

**対処**: `L = M` の場合を扱う補題は、**一般の `{L M} [Algebra L M]` の側に置く**。

```lean
theorem asIdeal_resHOS (V) :
    (resHOS (L := L) V).asIdeal = V.asIdeal.comap (algebraMap (𝓞 L) (𝓞 M)) := rfl
```

そこでは `((algebraMap (𝓞 L) (𝓞 M)) x : M) = algebraMap L M (x : L)` が `rfl` である
（mathlib の `RingOfIntegers.instAlgebra` はそう作られている）。
★同じ式を `L = M` に代入した文脈で**書き直す**と `Algebra.id` を拾って壊れる。

## `bernoulli` の値は `decide` でも `norm_num [bernoulli]` でも出ない

**失敗形**:

```lean
example : (bernoulli 4 : ℚ) = -1/30 := by decide
-- Decidable インスタンスが `bernoulli' 4` の再帰で止まる
example : (bernoulli 6 : ℚ) = 1/42 := by norm_num [bernoulli]
-- `bernoulli.eq_1` と `bernoulli'_def` が looping simp theorem になり maxRecDepth
```

★原因: `bernoulli'` は `bernoulli' n = 1 - ∑_{k<n} C(n,k)/(n-k+1) * bernoulli' k` という
**自分自身を含む** well-founded 再帰。simp に渡すと展開が止まらない。

**対処**: `bernoulli'_def` を **1 回だけ `rw`** して、下の値は**既に証明した補題**として
`norm_num` に渡す。`Nat.choose` も明示的に渡さないと `Nat.choose 4 2` が残る。

```lean
theorem bernoulli'_four : (bernoulli' 4 : ℚ) = -1/30 := by
  rw [bernoulli'_def]
  norm_num [Finset.sum_range_succ, bernoulli'_zero, bernoulli'_one, bernoulli'_two,
    bernoulli'_three, Nat.choose]
```

`n = 3, 4, 5, 6` と 1 段ずつ積み上げる。`bernoulli n` へは
`bernoulli_eq_bernoulli'_of_ne_one`(`n ≠ 1`)で移る。
★`riemannZeta 6` も mathlib に無いので、`riemannZeta_two_mul_nat` と `B₆` から作る
(`riemannZeta_two`・`riemannZeta_four` はある)。

## `Ring.inverse` は環準同型と交換しない(単元でない限り)

**失敗形**: `Ideal.Quotient.mk` や特殊化 `φ : R →+* R'` を `tateXterm t = t * Ring.inverse (1-t)^2`
のような式に当てて、`φ (Ring.inverse x) = Ring.inverse (φ x)` を暗黙に使ってしまう。

★`Ring.inverse x` は `x` が単元でないとき**既定値 `0`** を返す。`φ x` が単元になっても
`φ 0 = 0 ≠ Ring.inverse (φ x)` なので**一般には成り立たない**。

**対処**: `IsUnit x` を仮定に持つ。

```lean
theorem map_ring_inverse (φ : R →+* R') {x : R} (hx : IsUnit x) :
    φ (Ring.inverse x) = Ring.inverse (φ x) := by
  have hux : IsUnit (φ x) := hx.map φ
  have h : φ (Ring.inverse x) * φ x = 1 := by
    rw [← map_mul, Ring.inverse_mul_cancel _ hx, map_one]
  calc φ (Ring.inverse x) = φ (Ring.inverse x) * (φ x * Ring.inverse (φ x)) := by
        rw [Ring.mul_inverse_cancel _ hux, mul_one]
    _ = (φ (Ring.inverse x) * φ x) * Ring.inverse (φ x) := by ring
    _ = Ring.inverse (φ x) := by rw [h, one_mul]
```

★`Ring.eq_inverse_of_mul_eq_one_left` は**無い**ので上の `calc` で作る。
★★これが理由で、`Ring.inverse` を含む式の mod `I` 議論では `Ideal.Quotient` を使わず
**差の分解**(`Y²−g² = (Y−g)(Y+g)` 等)で処理する。

## `MvPolynomial.X` の素元性は mathlib に無い / `IsCoprime` と `IsRelPrime` を混同しない

**失敗形 1**: `Prime (MvPolynomial.X 0 : MvPolynomial (Fin 2) ℤ)` を `exact?` で探しても出ない
(`MvPolynomial.prime_X` も `MvPolynomial.irreducible_X` も存在しない)。

**対処**: `MvPolynomial.finSuccEquiv` で `Polynomial` に移す。

```lean
theorem prime_univA : Prime (X 0 : MvPolynomial (Fin 2) ℤ) := by
  have hp : Prime (Polynomial.X : Polynomial (MvPolynomial (Fin 1) ℤ)) := Polynomial.prime_X
  rw [← MvPolynomial.finSuccEquiv_X_zero (R := ℤ) (n := 1)] at hp
  exact (MulEquiv.prime_iff (M := MvPolynomial (Fin 2) ℤ)
    (MvPolynomial.finSuccEquiv ℤ 1 : MvPolynomial (Fin 2) ℤ ≃ₐ[ℤ] _).toRingEquiv).1 hp
```

★`X 1` は `finSuccEquiv_X_succ` で `Polynomial.C (X 0)` になり、`Polynomial.prime_C_iff` で落とす。
★★`MulEquiv.prime_iff` に `.toRingEquiv` を渡すときは `(M := …)` を明示しないと単一化に失敗する。

**失敗形 2**: `IsCoprime (X 0) (X 1)` を示そうとする。

★これは**成り立たない**——`IsCoprime` は Bezout の意味(`∃ a b, a*x + b*y = 1`)であり、
`(X 0, X 1)` は真のイデアルなので 1 を生成しない。

**対処**: `IsRelPrime`(共通の非単元因子が無い)を使う。`IsRelPrime.mul_dvd`
(`DecompositionMonoid` が要る。UFD なら在る)と `IsRelPrime.pow` で

    x^n ∣ P ∧ y^n ∣ P → (x*y)^n ∣ P

が出る。`IsRelPrime x y` 自体は「`x` が素元」＋「`x ∤ y`」から手で作る。

## 繊維で括り直すとき——`Equiv.sigmaFiberEquiv` に任せ、繊維ごとの同値だけ作る

**失敗形**: `ℕ × ℕ ≃ Σ N, ((N+1).divisorsAntidiagonal : Finset _)` を直接作ろうとする。
`right_inv` で `Sigma.mk.injEq` を書き換えると**第二成分が依存型**なので motive が壊れる。

**対処**: 繊維分解そのものは mathlib の `Equiv.sigmaFiberEquiv f` に任せ、
**繊維ごとの非依存な同値**だけを手で作る。

```lean
def fiberEquiv (N : ℕ) :
    {p : ℕ × ℕ // (p.1 + 1) * (p.2 + 1) - 1 = N}
      ≃ ((N + 1).divisorsAntidiagonal : Finset (ℕ × ℕ)) := …   -- Subtype.ext だけで済む
```

そのうえで `Summable.tsum_sigma` → `Finset.tsum_subtype` と流す。
★`(Equiv.sigmaFiberEquiv f) p = p.2.val` は `rfl`。
★★`Equiv.tsum_eq` は `∑' c, g (e c) = ∑' b, g b`。`rw [← e.tsum_eq g]` のあと
`exact tsum_congr fun c => rfl` で beta 差を潰す。

## `omega` は `(a, b).1 * (a, b).2` と `a * b` を別の原子として扱う

`show a * b = …` で形を揃えてから `omega` に渡すこと。
同様に `Nat.sum_div_divisors n id` は `∑ d, id (n / d)` の形なので、
`∑ d, n / d` とは構文的に一致しない——`show … = ∑ d, id (n / d) from rfl` を挟む。

## 三角不等式を項式で繋ぐとメタ変数が決まらない

**失敗形**:

```lean
have hbound := (norm_add_le _ _).trans (add_le_add (…) (le_of_eq (norm_mul _ _)))
-- don't know how to synthesize implicit argument `a` / `b` / `c` …
```

★`norm_mul _ _` の `_` は上流の `_` からは決まらない。`have` に型注釈が無いと
全部メタ変数のまま残る。

**対処**: 組み立てを**独立した補題**として型を書き切る。

```lean
theorem norm_three_comb (Au Aw As : ℂ) :
    ‖Au + Aw + (-2 : ℂ) * As‖ ≤ ‖Au‖ + ‖Aw‖ + 2 * ‖As‖ := by
  have h1 := norm_add_le (Au + Aw) ((-2 : ℂ) * As)
  have h2 := norm_add_le Au Aw
  have h3 : ‖(-2 : ℂ) * As‖ = 2 * ‖As‖ := by rw [norm_mul]; simp
  linarith
```

そのうえで本体は `have hcomb := norm_three_comb <式> <式> <式>` と当て、最後に `linarith`。
★★`linarith` に渡す形にしておくと、係数の帳尻(`8‖q‖^{n+2} ≤ (4‖q‖)^{n+1}` など)も
同じ `linarith` で片付く。

## 帰納法で回る形にするため、仮定はあえて弱くする(`w ≠ 0` に限る)

`‖f(w)‖ ≤ C‖w‖^m`(小さい `w`)から `X^m ∣ f` を出す帰納法で、
仮定を `∀ w, ‖w‖ < r → …`(`w = 0` を含む)と書くと**楽に見えるが帰納法が回らない**。

★`f = X * g` と割ったあと `g` について同じ形の仮定が要るが、`w = 0` での評価は
`f` の仮定からは出ない(`‖g(0)‖ ≤ C·0^m` は `‖f(0)‖ ≤ C·0^{m+1}` から導けない)。

**対処**: 仮定を `∀ w, w ≠ 0 → ‖w‖ < r → …` と**弱めて**おく。すると `f` でも `g` でも
同じ形が保たれ、`w = 0` の評価は毎段「連続性 + `𝓝[≠] 0` の極限」で作る。

```lean
theorem eval_zero_le_of_bound (g : Polynomial ℂ) (C r : ℝ) (hr : 0 < r) (m : ℕ)
    (h : ∀ w : ℂ, w ≠ 0 → ‖w‖ < r → ‖g.eval w‖ ≤ C * ‖w‖ ^ m) :
    ‖g.eval 0‖ ≤ C * ‖(0 : ℂ)‖ ^ m := by
  …  -- le_of_tendsto_of_tendsto + g.continuous.continuousAt + self_mem_nhdsWithin
```

★`le_of_tendsto_of_tendsto` は `f ≤ᶠ[b] g`(eventually)を取る。
`filter_upwards [self_mem_nhdsWithin, hball.filter_mono nhdsWithin_le_nhds]` で供給する。

## 同じ集合を二度書くと `isDefEq` が爆発する

**失敗形**:

```lean
refine poly_eq_zero_of_infinite_zeros _
  (Set.range fun z : UpperHalfPlane => Complex.exp (2 * ↑π * I * (z : ℂ)))
  infinite_exp_range ?_
-- (deterministic) timeout at `isDefEq` (maxHeartbeats 2000000 でも落ちる)
```

★補題 `infinite_exp_range` の集合と、ここで書いた集合を単一化しようとして、
`UpperHalfPlane` の coercion の展開で爆発する。

**対処**: **集合を書かず、補題側から推論させる**。

```lean
refine poly_eq_zero_of_infinite_zeros _ _ infinite_exp_range ?_   -- 0.03 秒
```

★★同型の引数を「念のため明示する」のは、Lean では**逆効果になることがある**。
既に補題が持っている形は、補題に決めさせるのが速い。


## `AnalyticAt.comp` は合成先を勝手に別の形に分解する

**失敗形**:

```lean
theorem analyticAt_shiftP (L : PeriodPair) (w s : ℂ) (h : s + w ∉ L.lattice) :
    AnalyticAt ℂ (fun u => L.weierstrassP (u + w)) s :=
  (L.analyticOnNhd_weierstrassP (s + w) h).comp (analyticAt_id.add analyticAt_const)
-- Type mismatch: has type AnalyticAt ℂ (℘[L] ∘ HAdd.hAdd s) w
--                but is expected to have type AnalyticAt ℂ (shiftP L w) s
```

★`AnalyticAt.comp : AnalyticAt g (f x) → AnalyticAt f x → AnalyticAt (g ∘ f) x` に
`AnalyticAt ℘ (s + w)` を渡すと、エラボレータは `f x` を `HAdd.hAdd s w` と読んで
**`f := HAdd.hAdd s`、`x := w`** と分解してしまう(欲しいのは `f := (· + w)`、`x := s`)。
`s + w` は `f x` として 2 通りに読めるので、先に来た方が選ばれる。

**対処**: **`f` と `x` を名前つき引数で明示する**。

```lean
AnalyticAt.comp (f := fun u : ℂ => u + w) (x := s)
  (L.analyticOnNhd_weierstrassP (s + w) h) (analyticAt_id.add analyticAt_const)
```

★★「集合を明示すると爆発する」場合(前項)と逆で、**合成の分解は明示しないと決まらない**。
分かれ目は「補題側が既にその形を持っているか」。持っていれば任せ、二通りに読めるなら明示する。

## `deriv (shiftP L w)` は `rw` で開かない

`noncomputable def shiftP L w := fun s => L.weierstrassP (s + w)` に対して

```lean
rw [shiftP]        -- Failed to rewrite using equation theorems for `shiftP`
rw [deriv_shiftP]  -- 続く rw が (fun s => ...) (z - l₀) の形で止まる
```

**対処**: 定義を開くときは `show`、書き換えたあとにベータ簡約が要るときは `simp only`。

```lean
show deriv (fun u : ℂ => L.weierstrassP (u + w)) s = _   -- 定義を開く
simp only [deriv_shiftP]                                  -- 開いてベータ簡約まで
```


## 3 重の `Polynomial` は `eval₂Hom` の行き先が推論できない

**失敗形**:

```lean
noncomputable def toPP3 : CollBase →+* Polynomial (Polynomial (Polynomial ℤ)) :=
  MvPolynomial.eval₂Hom
    (((Polynomial.C : Polynomial (Polynomial ℤ) →+* Polynomial (Polynomial (Polynomial ℤ))).comp
      (Polynomial.C : Polynomial ℤ →+* Polynomial (Polynomial ℤ))).comp
      (Polynomial.C : ℤ →+* Polynomial ℤ)) ![...]
-- failed to synthesize instance of type class
--   CommSemiring (Polynomial (Polynomial (Polynomial ?m.53)))
```

★戻り値の型注釈は `def` の側にあるのに、`eval₂Hom` の `S₁` は先に決まらない。
2 重(`Polynomial (Polynomial ℤ)`)までは通るが、3 重で `?m` が残ってインスタンス
探索が落ちる。

**対処**: **`abbrev` で名前をつけて `S₁` を名前つき引数で渡す**。

```lean
abbrev PPP : Type := Polynomial (Polynomial (Polynomial ℤ))

noncomputable def toPP3 : CollBase →+* PPP :=
  MvPolynomial.eval₂Hom (S₁ := PPP) (...) ![...]
```

★同じ形は `Polynomial.eval₂RingHom` の入れ子にも起きうる。深い塔を作るときは
中間の型に名前をつけておくと、エラーも読めるようになる。

## `Kˣ ⧸ Subgroup.zpowers Q` で `rw [one_mul]` が当たらない

**失敗形**

```
h : 1 * c = 1        (c : Kˣ ⧸ Subgroup.zpowers Q)
rw [one_mul] at h
-- Did not find an occurrence of the pattern 1 * ?a in the target expression 1 * c = 1
```

`simp`・`simpa`・`group` も「made no progress」で止まる。★単位群の商では
`MulOneClass` の実例が 2 経路(`QuotientGroup` 由来と `CommGroup` 由来)で来るため、
`rw` の統一が構文的に失敗する。一般の `G ⧸ N` では起きない。

**直し方——項の水準で書く(defeq で通る)**

```
exact (one_mul c).symm.trans h      -- 1 * c = 1  から  c = 1
exact (mul_one c).symm.trans h      -- c * 1 = 1  から  c = 1
```

★★`rw`/`simp` は構文照合、`exact` は defeq 照合。**実例のダイヤモンドは後者なら抜ける**。

## `R` が結論にしか現れない補題は `(R := R)` を明示する

`theorem foo (W : WeierstrassCurve K) [IsIntegral R W] (h : …) : 0 ≤ vAdd (tateDvrVal R K) …`
のように `R` がインスタンス引数と結論だけに現れる補題を `foo (C • W) hΔ` と当てると
`typeclass instance problem is stuck / IsFractionRing ?m K` で止まる。
インスタンス探索が結論より先に走るため。`foo (R := R) (C • W) hΔ` と書けば通る。

## 在庫は「名前」でなく「概念」で引く

新しい定義や補題を書く前に、付けようとしている**名前**ではなく
**何を証明しようとしているか**で grep する。名前で引くと自分の命名だけを見て
既存の同内容を見落とす(2026-08-26: Ward の定理を第 58 で証明済みなのに
再導出してしまった。`Somos`/`Eds`/`normEDS` のどれかで引けば一発だった)。


### ★★在庫を引くときは「結論の数字」で打つ(2026-08-26、違反 2 回目)

失敗形: `E₄³ − E₆² = 1728Δ` を自分で組み立てようとし、
`sturm_bound_levelOne`・`qExpansion_mul` という**使うつもりの道具の名前**で grep した。
→ 道具は全部見つかったので「結論は無い」と思い込んだが、結論は
`ModularForms/LevelOne/GradedRing.lean` に `discriminant_eq_E₄_cube_sub_E₆_sq` としてあった。

直し方: **結論に現れるリテラル**(`1728`)で grep する。1 行で当たる。
一般に、探すべきは「道具の名前」ではなく「結論の形」である
——数字・係数・特徴的な記号は名前より強い手がかりである。


### ★★引用の直後の行に `=>` を置かない(2026-08-26)

失敗形: docstring を

    原文 (GenEll p.17):
    > Proposition 3.4. ... For any -/
    theorem foo : Tendsto (fun τ : ℍ => ‖jFun τ‖) atImInfty atTop := by

と書いたら `node tools/check.mjs` が
「逐語が GenEll 物理 p.17 に見つからない(61/160 文字まで一致)」で落ちた。

原因: `check.mjs` の `QUOTE_RE` は引用本文を
`((?:[^
]*>[^
]*
)+)`——つまり **`>` を含む行が続く限り**とっている。
★Lean の `fun x => …` は `>` を含むので、**引用の次の行が式ごと食い込まれる**。

直し方: 引用行を docstring の最後にしない
——引用のあとに**空行 + 1 行の地の文**を置いてから `-/` で閉じる。
★直後の行に `>` が無ければ(普通の `theorem foo : Bar := by`)問題は起きないので、
**`=>`・`≥`・`->` を含む宣言行の直前だけ**気をつければよい。


### ★Python でアストラル面の文字を書くときは 8 桁エスケープ(2026-08-26、2 回見た)

失敗形: 記録を書く Python の中で、基本領域の記号(U+1D49F)を
サロゲート対(`\ud835` + `\udc9f`)で書いたら
`UnicodeEncodeError: surrogates not allowed` で落ちた。同じ穴に 2 度落ちた。

直し方: U+10000 以上の文字は **`\U0001D49F` の 8 桁形式**で書く。
★Python 3 の str はサロゲート対を結合しないので、UTF-8 への書き出しで必ず失敗する。
★★別の道: その文字を変数に入れて連結する(`D = u"\U0001D49F"` を `+` で繋ぐ)。

★★★2026-08-28 追記: **8 桁形式でも「どの O か」を間違える**。
本プロジェクトの整数環は U+1D4DE(𝓞 = MATHEMATICAL BOLD SCRIPT CAPITAL O)であって
U+1D4AA(𝒪 = MATHEMATICAL SCRIPT CAPITAL O)ではない。
★見た目がほぼ同じなので、書いた直後の目視では気づかない。

★★★★★2026-09-06 追記(**3 回目。損害の本体が抜けていた**)。

★★**落ちるだけでは済まない——ファイルが 0 行になる。**
`io.open(path, "w")` は**開いた時点で切り詰める**ので、
`.write()` が `UnicodeEncodeError` で落ちると、中身が消えたファイルだけが残る。
実例: `Skeleton/PGC/Section1.lean` が 208 行 → **0 行**(`git checkout --` で復旧)。
過去にも `VeluSemistableJ.lean` が 117 行 → 0 行になっている。

★**安全な型**(符号化を先に済ませ、一時ファイル経由で差し替える):

```python
out = s.replace(old, new)
data = out.encode("utf-8")          # ★ここで落ちれば原本は無傷
fd, tmp = tempfile.mkstemp(dir=os.path.dirname(P), suffix=".tmp")
with os.fdopen(fd, "wb") as f:
    f.write(data)
os.replace(tmp, P)                  # ★差し替えは原子的
```

★この型の効き目は実証済み。この項目を追記するスクリプト自体が同じ罠を踏んだが、
先に `encode` していたので `lean-idioms.md`(5905 行)は**無傷だった**。

★**どちらの O かはファイルを測って決める**。見目では判別できないので、
書き込む先を `collections.Counter(ch for ch in s if ord(ch) > 0xFFFF)` で数える。
実測(2026-09-06): `Found/PGC/LubinTateReciprocityIsomorphism.lean` は
**U+1D4AA(𝒪) を 131 回**使っており、U+1D4DE は 0 回。ファイルごとに違う。
★★**エスケープを書かず、ファイルから既存の行をコピーする**のが確実。
やむを得ず書いたら書いた後に grep で誤った方の字を数え、0 を確認する
(ビルドは「Unknown identifier」で捕まえてくれるが、docstring に混ざると気づかない)。


## `Module.IsTorsionBySet.module` の SMul は diamond を作る(2026-08-28)

失敗形: `M` が `I` で消えることから `Module (R ⧸ I) M` を作り
(`htor.module`)、`Module.Finite.of_restrictScalars_finite` に渡そうとしたら

    failed to synthesize IsScalarTower R (R ⧸ I) M

文脈に `IsScalarTower` が**あるのに**拾わない。

原因: `Module.IsTorsionBySet` が作る `SMul` に 2 つの経路がある。

  ・`htor.hasSMul`                            ← `mk_smul` 補題が使う形
  ・`DistribMulAction.toDistribSMul.toSMul`    ← `htor.module` を instance にした後に
                                                 `•` が解決される形

★ defeq だが**形が違う**ので、`rw [htor.mk_smul]` が
「Did not find an occurrence」で落ちる。
同じ理由で `htor.isScalarTower` が作るインスタンスも
(`Submodule.Quotient.instSMul'` ベース)、`of_restrictScalars_finite` が要求する
(`instSMul` ベース)と別物になる。

★★見分け方: エラー本文に `@instHSMul … htor.hasSMul` と
`@instHSMul … DistribMulAction.toDistribSMul.toSMul` が**並んで出る**。

★★★**直し方(2026-08-28 解決): `haveI` ではなく `letI` を使う**。

    letI := htor.module
    letI := htor.isScalarTower (S := R)

`Module` インスタンスは **data を持つ**ので、`haveI` だと中身が失われる
(`haveI` は proof-irrelevant に扱う)。すると `isScalarTower` が
`htor.module` の中身を参照できず、別のインスタンスとして扱われる。

★★★★見つけ方: mathlib の `Algebra/Module/Torsion/Basic.lean:597` が
`(hM : Module.IsTorsionBySet R M I) : letI := hM.module` と書いている。
**定義側が `letI` を使っているなら、使う側も `letI`** である。

★★★★★一般化: **instance が data を持つなら `letI`**。
`Module` / `Algebra` / `SMul` は data、`Finite` / `IsDomain` は Prop なので `haveI` でよい。

## structure の中に `/-! … -/` を書くと、そこで structure が終わる

失敗形: 欄をグループ分けしようとして

    structure Foo where
      a : ℕ
      /-! ### ここからは §4 の分 -/
      b : ℕ

と書いたら `unexpected identifier; expected 'lemma'` が**ずっと先の行で**出た。

直し方: **`--` の行コメントにする**。`/-! -/` はモジュール/セクションのドキュメントであり、
宣言の中には置けない。★フィールドの `/-- -/` は問題ない。
★★エラー位置が離れるので、「さっき追加した `/-!` は無いか」を先に見ること。

## `Finset.sum_union_inter` に `(s := …)` で引数を渡せない

失敗形: `Finset.sum_union_inter (s := U) (t := V) (f := fun p => …)` が
`Invalid argument name \`s\`` で落ちる(引数名が `s✝` の形でしか無い)。

直し方: **型を書いて `:=` で受ける**。

    have hui : (∑ p ∈ U ∪ V, f p) + ∑ p ∈ U ∩ V, f p
        = (∑ p ∈ U, f p) + ∑ p ∈ V, f p := Finset.sum_union_inter

## `lake` は `lean/` の中でしか走らせない

失敗形: リポジトリ直下で `lake build --dir=lean` を走らせたところ、
別のトゥールチェインで全再ビルドが始まり、タイムアウトで殺された拍子に
mathlib の olean が 1 つ途中で切れて `failed to read file … incompatible header` になった。

直し方: **`lake exe cache get!`**(`lean/` の中で)。★olean を手で削除しないこと。
★★そもそも `cd /d/Math_ABC3/lean` してから `lake` を走らせる。`--dir` で代用しない。

## 包む定義が要るインスタンスは `haveI` ではなく `letI`

失敗形: `haveI : Algebra E L := (IntermediateField.inclusion h).toRingHom.toAlgebra` の後で
`IsScalarTower.of_algebraMap_eq (fun _ => rfl)` が

    (algebraMap ℚ L) q is not definitionally equal to
    (algebraMap E L) ((algebraMap ℚ E) q)

で落ちる。同じ理由で `Subtype.ext rfl` も落ちる。

直し方: **`letI` にする**。`haveI` は命題としてのみ保持するので定義的な中身が消える。
★mathlib 自身(`FieldTheory/IntermediateField/Algebraic.lean`)も `let _ :=` で置いている。
★★中間体の次数比較は自分で塔を組まず
`IntermediateField.finrank_le_of_le_right` / `finrank_le_of_le_left` を引く。

## mathlib のモジュール名は動く——`import` は Glob で実在を確かめてから書く

失敗形(2026-08-27): 記憶にある名前で `import` を書いて、

    error: no such file or directory … Mathlib\Data\Complex\Module.lean
    error: bad import 'Mathlib.Data.Complex.Module'

実測で動いていたもの:

| 覚えていた名前 | 実在する名前 |
|---|---|
| `Mathlib.FieldTheory.Adjoin` | `Mathlib.FieldTheory.IntermediateField.Adjoin.Basic` |
| `Mathlib.Data.Complex.Module` | `Mathlib.LinearAlgebra.Complex.Module` |

直し方: **Glob(`.lake/packages/mathlib/Mathlib/**/<名前>*.lean`)で実在を見てから書く**。
★`abc3-lean` の `lean_start` は「読めた」としか言わないので、
**`unknown namespace` が出たら import が無いことを疑う**(名前の綴り違いではない)。

★★`Algebra ℚ ℂ` のインスタンスは `Mathlib.Algebra.Algebra.Rat` が要る。
`ℂ` を import しても付いてこない——`failed to synthesize Algebra ℚ ℂ` はこれ。

## REPL の `addToEnv` は「エラーが 1 つでもあれば」何も積まない

失敗形: 3 つの定理をまとめて `lean_check(addToEnv: true)` に渡し、3 つめだけ落ちたので
直して再送したら、1 つめ・2 つめも `Unknown identifier` になっていた。

直し方: **落ちたら、その塊全部を作り直して送る**。
★★塊が大きいなら、そもそもファイルに書いて `lake build` したほうが速い
(REPL は 1 定理ずつの探りに使う)。

## `Subring.closure` を台にした型を圏論の構造に載せると核が発散する

失敗形(2026-08-27、第 368-369 ブロックで 2 回落ちた):

1. `def D : ℕ ⥤ CommRingCat where obj n := CommRingCat.of (ratTower n); ...`
   と**素朴に `Functor` の構造体を書き**、`map_id` / `map_comp` を `ext; rfl` で埋める
   —— エラボレーションは通るが **核判定が 100 秒を超えて落ちない**。
2. `def towerMk (n) (q) (hq) : D.obj n := ⟨q, hq⟩` という**補助定義を挟む**
   —— これも核判定が 95 秒。

原因: `ratTower n = Subring.closure {...}` なので `CommRingCat.of (ratTower n)` の
`CommRing` インスタンスの経路が長く、核での defeq が重い。

直し方:

* `Functor.ofSequence` / `NatTrans.ofSequence`(`Mathlib/CategoryTheory/Functor/OfSequence.lean`)
  に寄せる —— `n ⟶ n+1` の射だけ与えれば `map_id` / `map_comp` は mathlib が持つ。**0.07 秒**。
* 補助定義を挟まず `show ratTower n from ⟨q, hq⟩` と**直に書く** —— **0.05 秒**。

★根本的に直すなら `closure` でなく**明示的な `carrier`** で部分環を定義する。
★★同じ穴: 「エラボレーションが速い」ことは「核が速い」ことを意味しない。
`lean_check` が 120 秒で背景に落ちたら、まずこの形を疑う。

## 関手の合成で作った対象は `rw`/`simp` を止める（`pullback` と `Over.forget`）

失敗形(2026-08-27、第 374-375 ブロック):

`D f := overRatTowerDiagram ⋙ Over.pullback f ⋙ Over.forget X` と定義すると
`(D f).obj i` は **`pullback (over.obj i).hom f` と defeq だが構文的に違う**。
その結果:

* `pullback.hom_ext` を当てると、生成される目標の項が
  「外側は関手合成の型・内側は生の `pullback` の型」という**混在**になる
* `Category.assoc` すら「pattern が見つからない」で落ちる
  （エラー末尾に `Note: The target expression is not type-correct under the
  instances transparency level` が出るのが目印）
* `simp only [..., pullback.lift_fst]` も同じ理由で発火しない

★型注釈を足して目標の**外側**を生の `pullback` に固定すると `Category.assoc` は通る:

    ((a ≫ b : (D f).obj m ⟶ pullback (over.obj i).hom f') = (c ≫ d : ...))

★★しかし**内側**（`pullback.lift` の引数の型）は依然として混在なので、
`pullback.lift_fst` はまだ発火しない。

直し方: **関手を合成で作らず、`obj` を生の形で直接書く**。

    noncomputable def D (f) : ℕᵒᵖ ⥤ Scheme where
      obj i := pullback (overRatTowerDiagram.obj i).hom f
      map h := pullback.map _ _ _ _ (over.map h).left (𝟙 _) (𝟙 _) _ _

★極限は合成版で取り、`Iso` で生版へ移す（`IsLimit.ofIsoLimit`）。
★★★教訓は「配管」の一般則と同じ: **defeq は `rw`/`simp` を助けない。**
構文をそろえるのは設計の仕事である。

## `Scheme.Spec.obj (op R)` と `Spec R` は別物として扱われる

★**2026-08-27、第 371-380 ブロックで 5 ブロックぶんの摩擦の正体がこれだった。**

mathlib の `AlgebraicGeometry.Spec (R : CommRingCat) : Scheme` は
`Scheme.Spec.obj (op R)` の `def` である。★**defeq だが構文的に別物**なので:

* `specZIsTerminal.from X : X ⟶ Spec (CommRingCat.of ℤ)`（`Spec` の形）
* 自分で `f : X ⟶ Scheme.Spec.obj (Opposite.op (CommRingCat.of ℤ))` と書くと（関手の形）

この 2 つを混ぜた瞬間、`rw`/`simp` が**すべて**落ちる。出るエラーは

    Note: The target expression is not type-correct under the `instances` transparency level

で、`Category.assoc` すら「pattern が見つからない」と言われる。
★★真の原因は最後の `Full error:` に出る:

    f has type X ⟶ Scheme.Spec.obj (Opposite.op (CommRingCat.of ℤ))
    but is expected to have type X ⟶ Spec (CommRingCat.of ℤ)

直し方: **`Spec R` の形に統一する**。`Scheme.Spec` を明示的に書くのは
関手性（`Scheme.Spec.map`、随伴、`mapCone`）が要るときだけにし、
対象は `Spec R` で書く。

★★★同じ穴の一般形は前節と同じ——**defeq は `rw`/`simp` を助けない**。
ただし今回のように**エラーの表面（`Category.assoc` が効かない）と原因（対象の書き方）が
遠い**ことがあるので、`Full error:` の末尾まで読むこと。

---

## `0_Source/*.txt` は逐語に使えない——`-enc UTF-8` の有無で記号が消える

2026-08-27、[Stacks] Lemma 29.40.4 で実際に読み違えた。

`0_Source/<論文>.txt` は `pdftotext` の**古い引き方**（`-enc UTF-8` なし）で作られており、
`→` `≥` `≫` `⊗` `≅` がすべて**空白に落ちている**。そのため

    (3) Ld is f -very ample for some d  1,
    (4) Ld is f -very ample for all  d  1,

の (3) が `d ≥ 1`、(4) が `d ≫ 1` であることが**テキスト層からは分からない**。
★「ある `d`」と「十分大きい全ての `d`」は主張の強さが違うので、これは実害である。

★★しかしこれは `pdftotext` の限界ではない。`tools/check.mjs` が実際に引く

    pdftotext -enc UTF-8 -f N -l N <pdf> -

は `→ ≥ ≫ ⊗ ≅` を**そのまま出力する**。落ちていたのは `.txt` の側だけだった。

**規約として引き出すこと:**

| 用途 | 使うもの |
|---|---|
| どのページに何があるか探す | `0_Source/*.txt`（速い。7654 ページを 1 秒で走査できる） |
| `.verbatim` / `原文 (...)` に写す | ★**必ず `pdftotext -enc UTF-8 -f N -l N` を引き直す** |
| 装飾（下線・書体・ハット）の確認 | ★★`pdftoppm -r 150` で目視。上の 2 つは代替しない |

★★★症状: ゲートが `S4 逐語が物理 p.N に見つからない(先頭 19/124 文字まで一致)` と言い、
`一致した末尾` が**ちょうど記号の直前で切れている**。
そこで止まったら、まず `.txt` から写していないかを疑うこと。

---

## `Scheme` を扱うファイルには `universe u` が要る

`variable {X Y S : Scheme.{u}}` と書くと

    error: unknown universe level `u`

`Scheme.{0}` で固定しないなら、`open` の後に `universe u` を 1 行入れる。
★mathlib の `AlgebraicGeometry` 側のファイルは全部そうしている。

---

## Python: アストラル面のエスケープと、台帳を壊す書き込み

2026-08-27 に **`ResearchPaper/mathlib-gap.json` を 2 回空にした**。

### 症状

    UnicodeEncodeError: 'utf-8' codec can't encode characters
    in position 15025-15026: surrogates not allowed

原因は `𝒪`（𝒪 を**サロゲート対を 2 つの `\u` で**書いたもの）。
Python 3 はこれを 1 文字に合成せず、**孤立サロゲート 2 個**として保持する。
`json.dumps(..., ensure_ascii=False)` は通るが、UTF-8 への符号化で落ちる。

★直し方: BMP 外は `\U0001D4AA` の 8 桁形で書く。`\uXXXX` を 2 つ並べない。

### ★★★本当に痛いのはここ

    io.open(p, 'w', encoding='utf-8').write(<エンコードで落ちる文字列>)

`open(p, 'w')` は**書く前にファイルを 0 バイトに切り詰める**。
`.write()` で例外が出た時点で、**元の内容は既に消えている**。

★★規約: 台帳を書き換えるときは**エンコードしてから開く**。

```python
out = (json.dumps(j, ensure_ascii=False, indent=1) + '\n').encode('utf-8')
io.open(p, 'wb').write(out)          # ここまで来れば必ず書ける
```

★★★落ちたら `git checkout -- <file>` で戻る。
**台帳を触る前にコミットしておくこと**が唯一の保険である。

### ★★★★2026-08-28 に **3 回目**をやった

規約を書いた翌日に `io.open(p,'w').write(json.dumps(...))` と書いて
また 0 バイトにした。★**「文字列を作る」と「ファイルを開く」を同じ行に書かない**。
コミット済みだったので `git checkout --` で戻った。

### ★★★★★2026-08-28 に **4 回目**をやった——今度は **Lean ファイル**

規約を自分で書いて、同じ日にまた 0 バイトにした。今回の違いは 2 つ:

1. ★壊したのは台帳ではなく **`Found/Arakelov/APicToSheaf.lean`**。
   規約を「台帳を書き換えるときは」と読んでしまい、Lean ファイルに適用しなかった。
   ★★**規約はすべてのファイル書き込みに適用する**。
2. ★`out = s` で**文字列は先に作った**のに `io.open(p,'w').write(out)` と書いた。
   ★★★**エンコードは `write` 時に起きる**ので、文字列を先に作っても意味がない。
   規約の正しい読み方は「**`.encode('utf-8')` を先に実行し、`'wb'` で開く**」である。

```python
data = s.encode('utf-8')     # ★ここで落ちてもファイルは無傷
f = io.open(p, 'wb'); f.write(data); f.close()
```

★★★★**もっと良い道: 非 BMP 文字をそもそも書かない**。
docstring や `.needs` の文字列なら `O_F` のように ASCII で代用できる。
Lean の識別子で必要なら**既存の行をコピーする**。

### ★docstring と宣言の間に `open ... in` を挟めない

失敗形: `/-- doc -/` の直後に `open scoped TensorProduct in` を置いたら

    unexpected token 'open'; expected 'lemma'

直し方: **`open ... in` を docstring の前に置く**。

```lean
open scoped TensorProduct in
/-- doc -/
theorem foo ...
```

### ★台帳を書き戻すときは `indent` を**元に合わせる**

`mathlib-gap.json` は **`indent=1`** である。`indent=2` で書き戻すと
全行が差分になり(864/859)、実際の変更 11 行が読めなくなる。
★書き戻したら必ず `git diff --numstat` を見て、行数が変更の大きさと合うか確かめる。

★★`len(out)` は**文字数**であってバイト数ではない。
日本語の台帳では `len` が減っても内容は減っていない(87515 バイト = 61495 文字)。

---

## 図式・錐が絡む `pullback.hom_ext` は「生の形」に**全部**揃える

2026-08-27、同型の spreading out（`Found/GenEll/IsoDescent.lean`）で 4 回落ちた。

`(baseChangeRatTowerDiagram f).obj n` と `pullback (overRatTowerDiagram.obj n).hom f` は
**defeq だが構文が違う**。`pullback.hom_ext` が作る目標では後者の形が現れるので、
前者の形の射を混ぜると `Category.assoc` すら発火せず、毎回

    Note: The target expression is not type-correct under the `instances` transparency level

が出る。`Full error:` を読むと `Quiver.Hom A ((D f').obj n)` と
`Quiver.Hom A (pullback ... f')` の食い違いだと分かる。

### ★★片方だけ揃えても駄目

最初は自作の射の**余域**だけを `pullback ...` にした。すると今度は
**域**（`(D f).map h'` の余域）が合わなくなり、同じエラーが場所を変えて出た。

### 直し方: 対象・射・錐の脚に別名を置いて 1 つの形に揃える

```lean
noncomputable abbrev bcObj (f) (n) : Scheme := pullback (overRatTowerDiagram.obj n).hom f
noncomputable abbrev bcPt  (f)     : Scheme := pullback (overRatTowerCone.pt).hom f
noncomputable def bcMap (f) (h : m ⟶ n) : bcObj f m ⟶ bcObj f n := (D f).map h
noncomputable def bcLeg (f) (n)         : bcPt f  ⟶ bcObj f n := (C f).π.app n
```

そのうえで `bcMap_fst` / `bcMap_snd` / `bcLeg_fst` / `bcLeg_snd` を `@[simp]` で置く。
★`bcLeg_*` の証明は `show (C f).π.app n ≫ _ = _` で**関手の形へ戻してから** `simp only` する
——補題の中では戻し、外では生、と役割を分ける。

★★★これで `rw [Category.assoc, bcDescHom_fst, …]` が普通に通るようになる。
**別名は中身を変えるためではなく、構文を 1 つに揃えて `rw` を通すために置く。**

### 付随して覚えたこと

* `reassoc_of% h` は `∀ {Z} (k : … ⟶ Z), …` を返す。**引数を 1 つ与える**必要がある:
  `(reassoc_of% keyAB) (pullback.snd _ f)`。
* `ℕᵒᵖ` は前順序の圏なので `Subsingleton (m ⟶ i)` が通る。
  余フィルターの `min` で 2 回落とした後、**2 本の道は自動的に等しい**
  （`Subsingleton.elim` → `rw` で片方に寄せる）。
* `have hki : IsCofiltered.min i j ⟶ i := IsCofiltered.minToLeft i j` のように
  `have` は Type 値にも使える。`obtain ⟨k, hki, hkj⟩ : ∃ … ` は
  `Prop` に潰れるので**射の取り出しには使えない**。

---

## `Γ(X, U)` の上では `ConcreteCategory.comp_apply` の `rw` が落ちる

2026-08-27、切断の降下（`Found/GenEll/SectionDescent.lean`）で踏んだ。

`Γ(X, U)` は `X.presheaf.obj (op U)` で、`X.presheaf` の型は
`TopCat.Presheaf CommRingCat X.toPresheafedSpace`
——これは `(Opens X)ᵒᵖ ⥤ CommRingCat` の **`def`** である。
そのため `rw [ConcreteCategory.comp_apply]` が

    Note: The target expression is not type-correct under the `instances` transparency level
    Full error: ... (D.obj l).presheaf has type TopCat.Presheaf CommRingCat ...
                but is expected to have type (Opens ...)ᵒᵖ ⥤ CommRingCat

で落ちる。`CommRingCat.hom_comp` / `RingHom.coe_comp` / `Function.comp_apply` は
**そもそも当たらない**（linter が「unused」と言う）。

★直し方: `rw` をやめて **`congrArg` で合成の外側を直接当てる**。

```lean
theorem app_map_comp_eq {D : ℕᵒᵖ ⥤ Scheme.{0}} {i j l : ℕᵒᵖ} (g : j ⟶ i) (k : l ⟶ j)
    {U : (D.obj i).Opens} {s t : Γ(D.obj i, U)}
    (h : (D.map g).app U s = (D.map g).app U t) :
    (D.map (k ≫ g)).app U s = (D.map (k ≫ g)).app U t := by
  rw [Functor.map_comp, Scheme.Hom.comp_app]
  exact congrArg (ConcreteCategory.hom (Scheme.Hom.app (D.map k) (D.map g ⁻¹ᵁ U))) h
```

`hom (A ≫ B) x` と `hom B (hom A x)` は**defeq なので `congrArg` は通る**。
`rw` はパターン照合なので通らない。★★「defeq は `rw` を助けないが `exact` は助ける」。

★★★もう 1 つ効いたこと: **図式 `D` を変数のままにする**。
`baseChangeRatTowerDiagram f` を直接書くと `D.obj l` が展開されて
同じエラーが別の場所で出る。**補題は一般の図式で書き、具体の図式は呼ぶ側で入れる。**


---

## `hz ▸` は「向きが合わない」と落ちる——`⊥` と `0` は同じ項ではない

2026-08-27、`Found/GenEll/VerticalBound.lean` で踏んだ。

```lean
have hJ : J ≠ 0 := fun hz => hI (le_bot_iff.mp (hz ▸ h))   -- ✕
```

    invalid `▸` notation, expected result type of cast is  I ≤ ⊥
    however, the equality hz of type J = 0 does not contain the expected result type

★`Ideal` では `0` と `⊥` が**同じ値だが同じ項ではない**ので、
`▸` も `rw [← hz]` もパターンが見つからない。

★★直し方: 書き換えをやめて**順序の推移で繋ぐ**。

```lean
have hJ : J ≠ 0 := fun hz => hI (le_bot_iff.mp (le_trans h hz.le))   -- ○
```

`hz.le : J ≤ 0` は `Eq.le` で取れ、`0` と `⊥` は順序の側では同一視される。

## `div_le_div_of_nonneg_right` は `0 ≤ c` を取る（`0 < c` ではない）

同じブロック。名前から `0 < c` を渡すと

    Application type mismatch: hd has type 0 < ↑(Module.finrank ℚ F)
    but is expected to have type 0 ≤ ↑(Module.finrank ℚ F)

★`c = 0` でも `a/0 = 0 ≤ 0 = b/0` で成立するので仮定が弱い。
★★`positivity` で `0 ≤ ↑(finrank …)` を出すのが一番短い。

## `Found/GenEll/` で `pullback.snd` を使うなら `open … Limits`

`open AlgebraicGeometry CategoryTheory` だけだと

    Unknown identifier `pullback.snd`

★`pullback` は `CategoryTheory.Limits` にある。ファイル先頭の `open` に
`Limits` を入れること——**型に `pullback` が出なくても、
`bcObj` のような `abbrev` を展開した先で必要になる**。

---

## 在庫確認は **新しい名前すべて**について引く

2026-08-27、`Found/GenEll/ComapMul.lean` を**二度書き**した。

平行セッションが 2026-08-17 に既に取っていた
`ideal_comap_eq_map_of_isAffine`（`ComapAffine.lean`）と
`comap_mul`（`ComapMul.lean`、**一般の射**）を知らずに、
`pullbackSpecIso` 経由で同じものを書き直し、
あまつさえ `Write` で相手のファイルを上書きした（`git checkout --` で復元）。

★原因は**在庫確認の部分実施**である。
新しい補題名を 5 つ確認したところで満足し、
後から足した `comap_mul_of_isAffine` / `ideal_comap_of_isAffine` を引かなかった。

★★手順（CLAUDE.md の「在庫」）:

```bash
node tools/decl-index.mjs                     # .cache/decl-index.txt を作る
for n in <新しい名前を全部>; do grep -c "\.$n\b" .cache/decl-index.txt; done
```

★★★**ファイル名でも引くこと**。`ComapMul.lean` が既にあることは
`ls lean/ABC3/Found/GenEll/ | grep -i comap` で 1 秒で分かった。

★★★★`Write` が「updated」と言ったら**既存ファイルである**。
「created」でなければ手を止めて `git log -- <path>` を見ること。

---

## `Point.map` は `rw` では当たらないが `refine … .trans` なら通る

2026-08-27、Tate 一意化の体拡大との両立（`Found/GaloisRep/TateCurveNatural.lean`）で踏んだ。

`tatePtPair` は `((tateCurveAt q hq).map (algebraMap R K)).toAffine.Point` の上にあり、
mathlib の `WeierstrassCurve.Affine.Point.map` は `(W'.baseChange F).Point` の上にある。
★`baseChange F = map (algebraMap R F)` は **`rfl`** だが、

```lean
rw [Point.map_some]     -- ✕
simp only [Point.map_some]  -- ✕（simp made no progress）
```

    Application type mismatch: … has type @Point K … ((tateCurveAt q hq).map (algebraMap R K)).toAffine
    but is expected to have type @Point K … (tateCurveAt q hq)⁄K

★★直し方: **`refine`/`exact` に載せる**。

```lean
refine (Point.map_some (S := R) φ (nonsingular_tateK a w q hq haw hwu hne hΔ)).trans ?_  -- ○
unfold tatePtPair
congr 1
```

★★★理由: **`rw` はパターン照合（`instances` 透明度）、`refine`/`exact` は
既定透明度の `isDefEq`** である。定義上等しいだけの型は `rw` を助けないが
`exact` は助ける——`Γ(X,U)` の項（本ファイル上方）と**同じ型の失敗形**である。

★「defeq なのに `rw` が落ちた」と思ったら、まず `refine (…).trans ?_` を試すこと。

## `unfold` の後はインスタンスが合流しない —— `show` で形を合わせる

**失敗形**（2026-08-27、`ProjectiveSpace.lean`）:

```lean
instance : IsProper (projSpaceOverSpec n R) := by
  unfold projSpaceOverSpec
  exact IsProper.mk
-- failed to synthesize IsSeparated (Proj.toSpecZero … ≫ Spec.map …)
```

ところが**同じ式を直に書けば通る**:

```lean
example : IsSeparated (Proj.toSpecZero … ≫ Spec.map …) := by infer_instance  -- OK
```

★`unfold` が展開した後の項は、インスタンス探索から見ると元の項と別物になる
（`letI` が残る・簡約段階が違う）。★★**直す形は `show`**:

```lean
instance : IsProper (projSpaceOverSpec n R) := by
  show IsProper (Proj.toSpecZero (MvPolynomial.homogeneousSubmodule (Fin (n + 1)) R)
      ≫ Spec.map (CommRingCat.ofHom (gradeZeroEquiv n R).toRingHom))
  exact IsProper.mk
```

★★★同じ理由で、定義の中で `letI : GradedAlgebra … := MvPolynomial.gradedAlgebra` と
書くより、ファイル冒頭で `attribute [local instance] MvPolynomial.gradedAlgebra` と
宣言するほうがよい——`letI` は展開後に項として残ってインスタンス探索を邪魔する。

## `IsProper` は合成の instance を持たない —— `IsProper.mk` を明示する

`IsSeparated`・`UniversallyClosed`・`LocallyOfFiniteType` はそれぞれ合成の instance を
持つが、**`IsProper (f ≫ g)` の instance は無い**（2026-08-27 実測）。
★`exact IsProper.mk` と書けば 3 つの親を instance 探索させられる。

## モノイダル関手の `ε`/`η` が `Iso.refl` なのに `rfl` が通らない —— `respectTransparency false`

**失敗形**（2026-08-27、`AmpleDef.lean`）:

```lean
example (s : (𝟙_ X.PresheafOfModules).obj (op (⊤ : X.Opens))) :
    trivValue (𝟙_ X.PresheafOfModules) ⊤ restrictPresheafUnit.symm s = s := rfl
-- Type mismatch: ?m = ?m vs …
```

`simp [trivValue, secOn, restrictPresheafUnit, Functor.Monoidal.εIso]` まで進めると
残る目標は `(Functor.OplaxMonoidal.η F).app o s = s` である。
★mathlib の `PresheafOfModules.pushforward₀OfCommRingCat` は
`εIso := Iso.refl _` / `μIso _ _ := Iso.refl _` で定義されているので**本当に恒等**だが、
既定の透明度では `rfl` が届かない。

★★**直す形**——mathlib 自身がその instance に付けているのと同じ option を付ける:

```lean
set_option backward.isDefEq.respectTransparency false in
theorem trivValue_unit_top (s : …) : trivValue … s = s := by
  simp [trivValue, secOn, restrictPresheafUnit, Functor.Monoidal.εIso]
  rfl
```

★★★一般に「mathlib 側の定義に `set_option backward.…` が付いていたら、
それを消費する側にも要る」と思ってよい。

## 前層加群の射の自然性は `PresheafOfModules.naturality_apply`

`e.hom.naturality` ではなく **`PresheafOfModules.naturality_apply`**（元の水準）である:

```lean
PresheafOfModules.naturality_apply (f : M₁ ⟶ M₂) (g : X ⟶ Y) (x : M₁.obj X) :
  f.app Y (M₁.map g x) = M₂.map g (f.app X x)
```

★`Over U` の site では `g` は `(Over.homMk (homOfLE h) : Over.mk (homOfLE h) ⟶ Over.mk (𝟙 U)).op`
と書く。★★`X` / `Y` は名前付き引数で明示しないと合わないことが多い。
★★★`exact hnat` か `exact hnat.symm` かはエラーメッセージの左右を見て決める
——`Eq.symm hnat` の表示が期待と**同じ向き**なら `exact hnat` が正しい。

## `Scheme.Modules.smul_Spec_def` は `rfl` だが `rw` も `#synth` も通らない

**測定**（2026-08-28、`Definition 1.1` の最後の 1 本を追ったとき）:

mathlib の `Mathlib/AlgebraicGeometry/Modules/Tilde.lean` には

```lean
instance : Module R Γ(M, U) :=
  inferInstanceAs <| Module R ((modulesSpecToSheaf.obj M).obj.obj (.op U))

lemma smul_Spec_def (r : R) (x : Γ(M, U)) :
    r • x = ((Spec R).presheaf.map U.leTop.op) ((Scheme.ΓSpecIso R).inv r) • x := rfl
```

があり、これが「`Γ(Spec ℂ, ⊤)` の `ℂ`-加群構造は `ΓSpecIso` 経由である」を与える。

★**ところが具体化すると instance が出ない**:

```lean
#synth Module ((CommRingCat.of ℂ) : Type)
  (Γ(unitModules (Spec (CommRingCat.of ℂ)), (⊤ : (Spec (CommRingCat.of ℂ)).Opens)) : Type)
-- failed to synthesize
```

★★したがって `rw [Scheme.Modules.smul_Spec_def]` も通らない
（`HSMul ↑(CommRingCat.of ℂ) ↑Γ(…) ?m` が出ない）。
`c : ℂ` を `c : ((CommRingCat.of ℂ) : Type)` に書き換えても同じ。

★★★**回避の方向**（未検証）: `moduleSpecΓFunctor` 側の項
（`(moduleSpecΓFunctor.obj M : Type)`）で書くと `•` は出る——
`arcFiber` はそちらの綴りなので、`ArcFiber.lean` の側から降りるほうが近い。

---

## `rw` は「同じ項の別の綴り」を見ない —— `have` で綴りを固定してから書き換える

2026-08-28、算術直線束の等長同型（`Found/Arakelov/AMetricIso.lean`）で 4 回落ちた。

### 症状

    Tactic `rewrite` failed: Did not find an occurrence of the pattern
      transUnit ?m.375 ?V (pullTriv ?φ ?V ?e) (pullTriv ?φ ?V ?e')
    in the target expression
      … transUnit (L * M).sheaf c.W (pullTriv (φ ⊗ᵢ ψ) c.W …) (pullTriv (φ ⊗ᵢ ψ) c.W …) …

★**パターンは目で見て一致している。** 落ちる理由は、補題の暗黙引数
`?L` が `?φ : ?L ≅ ?M` から決まり、`φ ⊗ᵢ ψ` の域は `L.sheaf ⊗ M.sheaf` と綴られるのに、
ゴールの項は `(L * M).sheaf` と綴られているからである。
★★`(L * M).sheaf = L.sheaf ⊗ M.sheaf` は **`rfl`** だが、`rw` は構文で照合する。

同じ落ち方が `metric` の欄でも出る:

    Application type mismatch: The argument tensorTriv … has type
      … ((restrictPresheafFunctor X c.W).obj (L.sheaf ⊗ M.sheaf)) …
    but is expected to have type
      … ((restrictPresheafFunctor X c.W).obj (L * M).sheaf) …

### ★★★直し方: `have` で**ゴールの綴り**に固定してから `rw`

```lean
have e1 : transUnit (L * M).sheaf c.W (pullTriv (φ ⊗ᵢ ψ) c.W (tensorTriv c.eA c.eB))
      (pullTriv (φ ⊗ᵢ ψ) c.W g)
    = transUnit (L' * M').sheaf c.W (tensorTriv c.eA c.eB) g :=
  transUnit_pullTriv (φ ⊗ᵢ ψ) c.W (tensorTriv c.eA c.eB) _   -- ← ここは defeq で通る
rw [e1] at e2
```

`have` の型注釈は `exact` と同じく**定義的等価まで**見るので、
補題の綴りとゴールの綴りの橋渡しはここで済む。

★最後の `rw [h1, h2, …]` も同じ理由で落ちることがある。
そのときは `Eq.trans` の連鎖（`exact h1.trans (h3.trans (hkey.trans h2.symm))`）に替える。
★★項の同一性だけを使う書き換えは `congrArg (fun t => f … t …) lemma` で作れる。

---

## 在庫確認を飛ばすと**同名で衝突してビルドが落ちる**（2026-08-28、1 日に 2 回）

`Found/` は 1 つの名前空間 `ABC3.Found.Arakelov` に 300 ファイル以上が入っている。
★新しい補題を書く前に**必ず**在庫を引くこと（CLAUDE.md 在庫）。

    node tools/decl-index.mjs        # .cache/decl-index.txt を作る
    grep -n "^\(theorem\|def\) *<名前>" .cache/decl-index.txt

★★`decl-index.txt` は**作った時点のスナップショット**なので、
同じセッションで足した宣言は載らない。**木も引くこと**:

    grep -rn "theorem <名前>\|def <名前>" lean/ABC3/

### 落ち方

1 回目（`evalOn_one`）——同じ import 木に無かったので**単体ビルドは通り**、
`Found.lean` に足した時点で衝突が出た。

2 回目（`arithGamma`）——`lake build <その 1 ファイル>` は通ったが、

    error: import ABC3.Found.Arakelov.Definition11 failed,
      environment already contains 'ABC3.Found.Arakelov.arithGamma.src'
      from ABC3.Found.Arakelov.AMetricHom

が `Found.lean` で出た。★**`.src` の名前も衝突する**——本体の名前だけ見ても足りない。

★★★**単体ビルドが通っただけで commit しないこと**。
`Found.lean` に import を足して `lake build`（全体）まで通してから commit する。

### ★★★★★3 回目は**台帳に嘘を書きかけた**（2026-08-28）

`PresheafOfModules.pullback` の `Monoidal` インスタンスが mathlib に無いのを見て、
台帳に `arakelov-pullback-monoidal` を新設した。★**在庫を引いていなかった**。

実際には `Found/Arakelov/PicSchemeDelta.lean` に

    instance pullbackPreOplax : (pullbackPre f).OplaxMonoidal

が `sorry` 無しで在り、さらに `isLocallyTrivial_pullbackPre`（`PicLTPull.lean`）の
証明の中に**自明化の輸送が明示的に書かれていた**（`bcIso`・`pullbackOnUnitIso`）。

★★同じ回に `schemeRingHom` / `schemePullback` / `pullbackFreeYonedaIso` を
書き下ろしたが、これらも `pullbackPhi` / `pullbackPre` / `pullbackFreeYonedaIso`
として**すべて在庫**であった（`pullbackFreeYonedaIso` は**名前まで同じ**）。

★★★**教訓**: 「mathlib に無い」を測ったら、**必ず続けて木を引く**。

    grep -rn "def <名前>\|theorem <名前>" lean/ABC3/
    grep -rn "<概念の日本語>" lean/ABC3/Found/*/*.lean | head

★★★★台帳に gap を新設する前は**とくに**引くこと——
台帳の嘘は Lean のビルドが捕まえてくれない。

### ★4 回目（2026-09-08、pGC `reciprocityUnits` の α-同変性）—— **`structure` だと、エラーが名指すのは自分が書いた名前ではない**

```
error: ABC3/Found.lean:1:0: import ABC3.Found.PGC.ReciprocityAlphaTransport failed, environment already contains 'ABC3.Found.PGC.ArtinDatum.mk.noConfusion' from ABC3.Found.PGC.Section3RealParameters
```

★書いたのは `structure ArtinDatum` だが、エラーが名指すのは
**自動生成される `ArtinDatum.mk.noConfusion`** である。
⇒ ★**このエラー文を `grep` するときは `.mk.noConfusion` を落として本体名で引くこと。**

★このときも `node tools/leanfile.mjs <その 1 ファイル>` と
`node tools/build.mjs ABC3.Found.PGC.ReciprocityAlphaTransport` は**どちらも通った**
（相手の `Section3RealParameters` を import していないため）。
★★**`lake build ABC3.Found` まで通すこと**——上の「落ち方」1・2 回目と同じ形が 3 度目である。

## `def` はインスタンス探索を塞ぐ（`abbrev` の後付けはできない、2026-08-28）

`gammaModPre R L := (ModuleCat.restrictScalars ρ).obj (L.obj (op ⊤))` と `def` で置くと、

    Module ↑Γ(Spec R, ⊤) ↑(gammaModPre R L)

が **見つからない**。書き下した形

    Module ↑Γ(Spec R, ⊤) ↑((ModuleCat.restrictScalars ρ).obj (L.obj (op ⊤)))

なら `ModuleCat.instModuleCarrierObjRestrictScalars` で **見つかる**。
★インスタンス探索は `def` の中身を見ない（`instances` 透明度）。

★★**後付けの `attribute [local reducible] foo` は拒否される**：

    failed to set `[local reducible]` for `foo`, recall that `[reducible]` affects
    the term indexing datastructures used by `simp` and type class resolution

★★★直し方は 2 つ：
1. 最初から `abbrev`（＝ `@[reducible] def`）で置く。
2. **書き下した形で定理を証明し、`def` 版は `rfl` 相当で受け直す**
   （`invertible_gammaRestrict` → `invertible_gammaModPre`、§9-788）。

★同じ穴が「戻り値の型」でも起きる：`pullSec f L ⊤` の行き先は
`op ((Opens.map f.base).obj ⊤)` の成分で、`op ⊤` とは `rfl` だがインスタンスは見つからない。
→ `pullSecTop` / `psiU` のように **戻り値・引数の型を `op ⊤` で宣言し直す**（§9-786、§9-787）。
`gammaSheafifyM`（§9-780）が最初の例である。

## 配管——`node -e "…"` の中のバッククォートは bash に食われる（2026-08-28）

台帳（`ResearchPaper/*.json`）を `node -e "…"` で書き換えたとき、
JS 文字列の中の **バッククォート付き Markdown**（``` `X` ```）が
**bash のコマンド置換**として実行され、その部分が**空文字に置き換わって**書き込まれた。

    $ node -e "… what: '★`X` が `ℤ`-固有 …' …"
    bash: X: command not found
    → 書き込まれたのは 「★ が -固有…」

★`node tools/check.mjs` は PASS する（Lean も JSON も壊れていない）ので、
**ゲートでは捕まらない**。書いた内容を読み返して初めて分かる。

★★直し方: **スクリプトを `.mjs` ファイルに書いて `node file.mjs` で走らせる**
（`$CLAUDE_JOB_DIR/tmp` に置く）。ヒアドキュメント（`<<'EOF'`、クォート付き）でもよい。
★★★同じ理由で `git commit -m "…`X`…"` も禁物である——`-F -` とヒアドキュメントを使う。

## 配管——Serre への道で踏んだ穴（2026-08-28、§9-808〜822）

### 1. `HomogeneousLocalization.NumDenSameDeg.num_add` / `deg_add` は `x` が明示引数

    num_add c1 c2   -- ✗ c1 が x(分母の submonoid)に食われる
    num_add (Submonoid.powers f) c1 c2   -- ○

★`{𝒜}` は暗黙、`x` は**明示**である。

### 2. `(c1 + c2).deg` は `num` の**型に現れる**

`rw [hdeg]` を先にすると motive が通らない。★**`num` を先に書き換える**こと。
同じ形は `exists_pow_of_numDenSameDeg`（`c.deg = k` を `c.num` の前で `rw` できない）でも出た。

### 3. `x_i^k ≠ 0` に `pow_ne_zero` は使えない

`R` 一般では `MvPolynomial` に零因子がありうるので `IsReduced` を要求される。
★`rw [MvPolynomial.X_pow_eq_monomial]; simp`（係数が `1` の単項式）で落ちる。

### 4. `def` で包んだスキームは `.Opens` が合わない

    (projSpace N R).Opens   -- ✗ Proj.basicOpen と構文的に合わない
    (Proj (…)).Opens        -- ○

★`projSpace` は `def`（semireducible）だからである。`Proj (…)` の綴りで書く。

### 5. 指数の等式は `▸` ではなく `subst`

`c.n ∣ L` から `L = c.n * k` を取り出したら **`subst`** する。
★`M^{⊗L}` と `M^{⊗(c.n·k)}` の型の食い違いはそれで消える（`▸` で運ぶ必要はない）。

### 6. `whiskerLeftIso` を使う（`◁` は射用）

    M ◁ (iso)              -- ✗ `◁` は射を取る
    whiskerLeftIso M (iso) -- ○

## 在庫——「無い」と書く前に `decl-index` を引く（2026-08-28）

`isLocallyTrivial_sheafify` の証明を読んで「層化した側の自明化には**名前が付いていない**」
と台帳に書いたが、**誤りだった**——`Found/Arakelov/SheafifyTriv.lean` に
`sheafifyTriv` / `sheafifyTrivOf` が既にあり、`transUnit_sheafifyTriv`・
`sheafifyTriv_restrict` まで揃っていた。

★原因は「証明の中を読んだだけで在庫を引かなかった」ことである。
★★CLAUDE.md の在庫の規律どおり、**まず `node tools/decl-index.mjs` を作って
`.cache/decl-index.txt` を grep する**。木を読むのは在庫を引いた後でよい。

## 同型の向き——`≅` は「どちらから」を必ず確かめる（2026-08-28）

    restrictPresheafUnit : 𝟙_ ≅ (restrictPresheafFunctor X U).obj (𝟙_)
                           ^^^^  ここが左

`(restrict).obj (𝟙_) ≅ 𝟙_` **ではない**。`tensorPowTriv` の基点を
`restrictPresheafUnit` と書いたため、`trivValue … = 1` が `rfl` でも `simp` でも
落ちず（残ゴールは `ε.app _ ((𝟙_).map _ 1) = 1` の形）、半日詰まった。
★正解は **`.symm`** を付けるだけであった。

★★合図は在庫にあった——`trivValue_unit_top`（`AmpleDef.lean`）が
`restrictPresheafUnit.symm` を使っている。**在庫の使用例を 1 つ読めば向きは分かる**。

★★★教訓を一般化すると: `≅` を引数に取る補題が `rfl` で落ちないときは、
**証明の中身を疑う前に同型の向きを `#check` する**。0.07 秒で分かる。

### ついでに——`V = ⊤` 限定の在庫は一般化できることが多い

`trivValue_unit_top` は `V = ⊤` 用だったが、証明（`simp [trivValue, secOn,
restrictPresheafUnit, Functor.Monoidal.εIso]; rfl`）は**そのまま任意の `V`** で通り、
右辺が `s` から `s|_V` に変わるだけであった（`trivValue_unit'`）。
★「⊤ でしか無いから使えない」と諦める前に、証明をコピーして一般形を試す。

## `⨆ U i = ⊤` を `rw` で入れると motive が壊れる（2026-08-28）

    rw [← hcov]   -- ✗ motive is not type correct

目標に `le_iSup U i : U i ≤ ⨆ U i` という**証明項が現れる**からである
（`homOfLE (le_iSup U i)` の引数）。`⨆ U i` を `⊤` に書き換えると、
その証明項の型が合わなくなる。

★**直し方は「同型で運ぶ」**である:

    have hup : (⊤ : X.Opens) ≤ ⨆ i, U i := le_of_eq hcov.symm
    -- 両向きの制限射 map (homOfLE hup).op / map (homOfLE le_top).op が互いに逆

★★`Opens` の射は `Subsingleton` なので、合成した射は `Subsingleton.elim` で
好きな射に取り替えられ、`A ⟶ A` は `𝟙 A` に潰れる（`map_id_apply`）。

### 付随——`PresheafOfModules.map` は `map_comp` が直接使えない

`restrictScalars` を挟むので `M.map g.op ≫ M.map f.op` が型検査を通らない。
★`show M.presheaf.map … ` で **`Ab` 値の前層に降りてから** `map_comp` を使う
（`Found/GenEll/SheafifyGlue.lean` の `map_map_apply`）。

## `haveI` で置いたインスタンスを解決が拾わないことがある（2026-08-28）

    haveI hspec : IsClosedImmersion (Spec.map …) := …
    exact IsClosedImmersion.comp _ _        -- ✗ failed to synthesize

`hspec` の型は目標と字面が一致しているのに instance 解決が拾わなかった。
★**明示的に渡す**と通る:

    exact @IsClosedImmersion.comp _ _ _ f g inferInstance hspec

★★`def` で包んだ射（ここでは `globalChartMorphism`）が絡むと起きやすい
——`show` で展開しても解決器は元の形を探しに行くことがある。
★★★「instance が見つからない」と言われたら、**まず `@` で手渡してみる**。

## 名前の衝突——プロジェクト側の `FinitePlace` が mathlib の同名を隠す（2026-08-28）

    namespace ABC3.Found.GenEll
    open NumberField
    …  ∏ᶠ v : FinitePlace K, ⨆ i, v (x i)   -- ✗ Function expected at v

`ABC3.Found.GenEll.FinitePlace ≝ IsDedekindDomain.HeightOneSpectrum (𝓞 F)`
（`Found/GenEll/ArithDiv.lean:65`）が **`NumberField.FinitePlace` を隠す**。
名前空間の中では自分の `abbrev` が優先されるからである。

★**直し方**: `NumberField.FinitePlace K` と**完全修飾で書く**。

★★合図は「`h` の側は `NumberField.FinitePlace` と表示されるのに、
目標の側は `FinitePlace` と表示される」ことである
——**表示が食い違ったら名前解決を疑う**。

★★★REPL では `open … in example` の形で書いていて通ったのに
ファイルでは通らなかった。名前空間の中かどうかで解決が変わるので、
**REPL の断片が通ってもファイルでもう一度ビルドする**。

## `|(n:ℝ)| = (n.natAbs : ℝ)` は `Int.cast_natAbs` では通らない（2026-08-28）

    rw [Int.cast_natAbs]        -- ✗ Unknown constant
    push_cast [Int.abs_eq_natAbs]  -- ✗ does nothing
    rw [← Int.cast_abs, Int.abs_eq_natAbs]; simp   -- ○

★先に `ℝ` の絶対値を `ℤ` の絶対値へ戻してから `Int.abs_eq_natAbs` を当てる。

## `aeval g (X j) = g j` は `rfl` ではない —— 依存型の輸送が要る（2026-08-28）

    hyperGen N R (Fin.succ i') = X i'                    -- ○ rfl
    hyperplaneHom N R (X (Fin.succ i')) = X i'           -- ✗ rfl でない

★`hyperplaneHom = aeval (hyperGen …)` で、`aeval_X` は `rfl` ではないからである
（`eval₂` は `Finsupp.sum` を通る）。

★★これが効くのは **`Away ℬ f` のように `f` が型に現れるとき**である:

    Away ℬ (hyperplaneHom N R (X i))   と   Away ℬ (X i')

は**定義的に等しくない**ので、`rw [hyperplaneHom_X_succ]` で
**型ごと書き換える**（motive は型検査を通る）か、`Eq.mpr` で運ぶ必要がある。

★★★添字の型にも注意: `awayEval N R i` は `i : Fin (N+1)` で
`MvPolynomial (Fin (N+1)) R` の上の写像である。
**終域側（`MvPolynomial (Fin N) R`）に使うには `N = M+1` の分解が要る。**

## `Away ℬ f` の `f` が型に現れる —— 変数に特殊化した補題は使えない

`HomogeneousLocalization.Away ℬ f` は **`f` を型の引数に持つ**。したがって

* `awayEval N R i`（`f = x_i` に特殊化）を作っておいても、
* 終域が `Away ℬ (hyperplaneHom (x_i))` の場面では**そのままでは使えない**

——`hyperplaneHom (x_{i+1}) = y_i` は**命題としては成り立つが `rfl` ではない**ので、
型が合わない。`rw`-in-type や `Eq.mpr` で運ぶのは苦しい。

**直し方**: 最初から `f` を**任意の次数 1 の斉次元**にしておく
（`awayCoordOf R f hf j`、`awayEvalOf R f hf`；`ABC3.Found.GenEll.AwayEvalGen`）。

**測定 (2026-08-28)**: 一般化しても証明は**一字も変わらなかった**
——`x_i` であることを使っていたのは「分母が単項式であること」ではなく
「**次数が 1 であること**」だけだったから。
`awayCoordOf R (X i) _ j = projCoord N R i j` は **`rfl`** なので、
既存の特殊形は一般形の定義的な場合として残る。

**系**: 次数付き環準同型 `g` についての四角
`Away.map g f ∘ awayEvalOf f = awayEvalOf (g f) ∘ g` に
「`g (C r) = C r`」の仮定は**要らない**——`awayEvalOf_mk` が**次数だけ**で
右辺を `Away.mk` に潰すからである。一般化した方が証明が短くなる例。

## `f.app U` は `f` に依存した型を持つ —— でも `ker` は持たない

`Scheme.Hom.app f U : Γ(Y,U) ⟶ Γ(X, f ⁻¹ᵁ U)` は**終域が `f` に依る**ので、
`h : f = g` があっても `f.app U = g.app U` とは**書けない**（型が違う）。
`congrArg (fun m => m.app U) h` は型エラーになる。

**直し方**: 核を取ってから比べる。`RingHom.ker (f.app U).hom : Ideal Γ(Y,U)` の型は
**`f` に依らない**ので、変数に対して `subst` できる:

```lean
theorem ker_app_congr {X Y : Scheme} {f g : X ⟶ Y} (h : f = g) (U : Y.Opens) :
    RingHom.ker (Scheme.Hom.app f U).hom = RingHom.ker (Scheme.Hom.app g U).hom := by
  subst h; rfl
```

同じ手が「開集合が型に現れる」場合にも効く——`chartA ⁻¹ᵁ (chartA ''ᵁ ⊤) = ⊤` は
`rw` では動かせないが、**開集合を変数として受け取り** `subst` すれば消える:

```lean
theorem ker_specMap_app_eqToHom (φ : B ⟶ B') (V : (Spec B).Opens) (h : (⊤ : _) = V) : … := by
  subst h; …
```

## `CommRingCat.hom_comp` の `rw` が「型が instances 透明度で正しくない」で通らない

`Γ(X, U)` を含む式で `rw [CommRingCat.hom_comp]` が

> Note: The target expression is not type-correct under the `instances` transparency level

で失敗することがある（`presheaf.obj` の型が `TopCat.Presheaf` と関手型の間で合わない）。

**直し方**: 宣言の直前に

```lean
set_option backward.isDefEq.respectTransparency false in
```

を置く。mathlib 自身が `Mathlib/AlgebraicGeometry/IdealSheaf/Basic.lean` の
`Hom.ker_apply` などで同じことをしている。
`simp only [CommRingCat.hom_comp]` でも駄目（「unused」と言われる）なので、
option を置くのが正解である。

## `rw [← f_hom]` が通らないのは結合の向きのせい（透明度の警告は副次的）

`a ≫ b ≫ c` は **`a ≫ (b ≫ c)`**（`≫` は右結合）なので、`a ≫ b` は**部分項ではない**。
`rw [← lemma]` で `a ≫ b` を畳もうとすると

> Did not find an occurrence of the pattern …
> Note: The target expression is not type-correct under the `instances` transparency level

が出る。**2 行目は副次的な症状**で、原因は 1 行目——結合の向きである。

**直し方**: 先に括り直す。

```lean
rw [← Category.assoc a, ← Scheme.Hom.appIso_hom f U]
```

`← Category.assoc` に**最初の射を明示的に渡す**と、どの箇所を括り直すかが決まる。
（2026-08-28、`Γ(Proj 𝒜, D₊(x_i)) ≅ A⁰_{x_i}` の打ち消しで実測）

## 因子の Green 関数は連続でない —— 仮定は「差の連続性」で置く

算術因子 `D̄ = (D, g)` の Green 関数 `g(p) = −log‖s(p)‖` は
**台 `|D|` の上で発散する**ので、`Continuous g` を仮定に置くと
**本物の因子には当たらない補題**ができる（空虚ではないが使えない）。

原典が言うのは「同じ直線束の**2 つの計量**の比が有界」であり、
そこで連続なのは**差** `D.green − E.green` の方である（特異性が打ち消し合う）。

**直し方**: 比較の補題は

```lean
(hcont : Continuous (fun p => D.green p - E.green p))
```

で書く（`ABC3.Found.GenEll.htArith_sub_abs_le_of_diff`）。
`D.divisor = E.divisor` なら有限側は打ち消し合い、`archADiv` は Green 関数について
線型なので、既存の一様評価がそのまま効く。

**測定 (2026-08-28)**: この形の「強すぎる仮定」は**証明が通るので気づきにくい**。
消費側で実際の対象（ここでは Fubini–Study 計量）を当てようとして初めて露見した。

## 在庫が「無い」と出たら、**主語を変えて引き直す**

`differentIdeal … = ⊤` の判定を探して

```
Algebra.isUnramified_iff_differentIdeal_eq_top   無い
differentIdeal_eq_top_iff                        無い
differentIdeal_self                              無い
```

で「mathlib に無い」と結論し、仮定として受けた（`§9-901`）。
**しかしそれは「その名前で無い」でしかなかった。**

`differentIdeal` を主語にした判定は確かに無いが、
**判別式を主語にすれば全部あった**（2026-08-28 実測）:

```lean
Algebra.finrank_eq_one_iff_bijective_algebraMap                    -- [L:K] = 1 ⟺ 全単射
NumberField.discr_eq_discr_of_algEquiv                             -- K ≃ₐ[ℚ] L ⟹ disc 一致
NumberField.natAbs_discr_eq_absNorm_differentIdeal_mul_natAbs_discr_pow
Ideal.absNorm_eq_one_iff                                           -- N(I) = 1 ⟺ I = ⊤
```

これで `[L:K] = 1 ⟹ 𝔡 = ⊤` が 10 行で出る（`differentIdeal_eq_top_of_finrank_eq_one`）。

**直し方**: 求めている命題を**同値な別の量で言い換えてから**引き直す。
`𝔡` なら `disc`・`absNorm`・`ramificationIdx`、
`高さ` なら `deg`・`absNorm`、`連結` なら `IsPreconnected`・`Irreducible` など。
`exact?` が timeout するのは「近い形も無い」ではなく「主語が違う」の合図でもある。

## 原文が「by working locally」と書いた段を大域へ持ち上げない

`[GenEll] Proposition 1.7` の elementary claim は
**`ℚ_p` の有限次拡大**についての主張であり、原文はその直前に
「Moreover, **by working locally**, we reduce immediately to …」と書いている。

これを数体（`NumberField.ringOfIntegers`）へ持ち上げると**偽**になる:

- `K = ℚ(ζ_3)`、`L = K(∛2)`、`p = 3`
- `2` の上の素点 `𝔮` は `e = 3`・剰余標数 `2`（`2 ∤ 3` なので**馴**）→ `v_𝔮(𝔡) = 2`
- `𝔮 ∤ 3` なので `v_𝔮(3^3) = 0` → `3^3 ∉ 𝔡`

局所体（剰余標数 `p`）なら `v(p) ≥ e(L/K)` なので馴の場合も `e−1 < v(p)` で通る。
**`p` と無関係な素点で馴分岐が起こり得るのが大域だけの現象**である。

**直し方**: 原文の「locally」「at a prime」「in a neighborhood」は
**形式化でも局所環・付値環のまま置く**。大域へ上げるのは、
上げても成り立つことを**素点ごとの不等式で確かめてから**にする。

**測定 (2026-08-28)**: 大域版は**条件付き定理としては通る**（仮定 `step` を持つ形）ので、
`lake build` では何も出ない。仮定を充足しようとして初めて露見した
——`§9-872`（Green 関数の連続性）と**同じ形の失敗**である。

## `set_option … in` は**ドキュメンテーション文字列より前**に置く

```lean
/-- doc -/
set_option maxHeartbeats 400000 in
theorem foo : True := trivial
```

は **パースエラー**である:

```
unexpected token 'set_option'; expected 'lemma'
```

正しい順序は逆:

```lean
set_option maxHeartbeats 400000 in
/-- doc -/
theorem foo : True := trivial
```

**測定 (2026-08-28)**: 同じファイルに 3 箇所同じ誤りがあったが、
`lake build` が報告したのは **3 番目だけ**だった
——パーサのエラー回復が先の 2 つを飲み込んでいる。
★「1 件だけ直せば通る」と思って直すと、次のビルドで次の 1 件が出る。
**同じ形が他にもないか、その場でまとめて grep すること。**
★2026-09-06 に 4 回目。`Skeleton/` の statement を**変えずに**
`linter.unusedSectionVars` を黙らせたくて、既にある docstring と `theorem` の
あいだに差し込んで踏んだ。★`omit [Inst] in` で消してはいけない
——型から instance を落とすので **statement の変更**になる。

## `f ∣_ U` を含むゴールで `rw [f = g]` は motive が型付かない

`morphismRestrict`（`f ∣_ U : f ⁻¹ᵁ U ⟶ U`）は**始域が `f` に依存する**ので、

```lean
have hfac : ψ' ≫ V.ι = ψ
rw [hfac]   -- motive is not type correct
```

は通らない（`fun _a => IsImmersion (_a ∣_ U)` の型が合わない）。

**直し方**: 汎化してから `rintro rfl`。

```lean
have hgen : ∀ χ, χ = ψ → IsImmersion (χ ∣_ U) := by
  rintro χ rfl
  …  -- ここでは χ が ψ に置き換わっている
exact hgen _ hfac
```

同じ形は `f ⁻¹ᵁ U`・`f.app U`・`Scheme.Opens.ι` を含むゴールでも起きる
（`§9-864` の `ker_app_congr` は同じ病に対する別の処方）。
（2026-08-28、`§9-916` で実測）

## 引用照合(`check.mjs`)—— 本文からコピーした逐語が落ちることがある

**失敗形**

```
NG  lean\ABC3\Found\GenEll\Thm21Chain.lean:84
    引用照合: 逐語が GenEll 物理 p.12 に見つからない(layout で 3/65 文字まで一致)
    次に来るはず: "factthat(i)=⇒(ii)isimmediatefromthedefinitions.Thus,itsuffic"
```

**原因** `.txt`(`ResearchPaper/0_Source/*.txt`)から目で写した逐語と、
`check.mjs` が `pdftotext -layout` で取り直した文字列が**別の字**になることがある。
数式記号（`⇒` `≲` `ϵ′` `ω`）や合字の周りで起きる。
★★★**訂正（2026-09-07、メタ第 21 回）**: 旧版の「次に来るはず」に出ていたのは
★**PDF 側ではなく我々の写しの方**だった。★★**本体はこれを PDF 側だと誤読して 4 往復失った。**

★現在の NG は **我々の写し** と **PDF の実物** を並べて出す（`我々の写し →` / `PDF の実物 →`）。
★★**変種は 3 つ**（`layout` / `default` / `raw`）あり、`★` が付いた変種が最も長く一致したもの。★**どれか 1 つに当たれば照合は通る**
（★実測: 逐語の **10%** は `raw` でしか当たっていない。★**3 変種は飾りではない**）。

★★★**引用を書く前に投影を見るのが最速**:
```
node tools/check.mjs --projection --paper <鍵> --find '<語>'
```
★探し語も `leanQuoteProjection` を通るので **docstring の形をそのまま貼ける**。
★実例（2026-09-07 の実測）: pGC は `ΓKab`（★下付き K が `ab` の**前**）/ `d=ef` / `→ Z → 0`（★ハット無し）/ `OK`、
Yoshida は `W(Kmf/K)`（★ハット無し）/ `[θ](j)`（★波括弧無し）/ `ρf,m`。

**直し方** その項目で**すでに通っている逐語**（多くは見出し行）に差し替える。
落とした逐語の内容は地の文に `` ` `` で括って残せばよい。

```lean
原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

★原文 p.12 は「`The fact that (i) ⟹ (ii) is immediate from the definitions.`」と書く。
```

★ページ番号は**その docstring の `原文 (Tag p.NN):` 行**から取られる（`.src` の
`pdfPage` ではない）。差し替えるときは両方を見ること。

**★コミット前に必ず `node tools/check.mjs` を見る。**
`lake build` が通っても引用照合は落ちる——2026-08-29 に NG 2 件のまま
コミットしてしまった（次のコミットで修復）。

## Python ヒアドキュメントの `\uXXXX` サロゲート —— **書きかけでファイルを壊す**

**失敗形**（2026-08-29、2 度目）

```
UnicodeEncodeError: 'utf-8' codec can't encode characters in position ...: surrogates not allowed
→ その直後の lake build で
error: ABC3/Skeleton/GenEll/Section3.lean:196:22: Unknown identifier `EllModuliData`
```

**原因** `𝔽`（`𝔽` = アストラル面 U+1D53D）のような**サロゲートペアの
エスケープ**を Python 文字列に書くと、`str` としては不正なまま作られ、
`io.open(...).write()` の時点で例外になる。
★**例外は書き込みの途中で起きうる**ので、**ファイルが壊れた状態で残る**。

**直し方**

* **`Edit` ツールを使う。** Lean ファイルの内容は Write/Edit で書く（CLAUDE.md の規則）。
  アストラル面の文字を Python 経由で流し込まない。
* どうしても Python を使うなら、`𝔽` ではなく**文字そのもの**（`𝔽`）を
  ソースに直接書く。エスケープに割らない。
* ★**壊れたら `git checkout -- <path>` で即座に戻す。** `git status` で
  他に触っていないことを確かめてから続ける。

★★CLAUDE.md の「配管」にも同じ形が登記されている（アストラル面エスケープ）。
**2 度目なのでここに失敗の"結末"（ファイルが壊れる）まで書いた。**

## 在庫の測り方 —— **grep だけの「無い」は弱い**

**失敗形**（2026-08-29）

`Mathlib/AlgebraicGeometry/` と `Mathlib/NumberTheory/` を grep して
「mathlib に Kummer 対応は無い」と `.needs` に書いた。
★**`Mathlib/FieldTheory/KummerExtension.lean` を見ていなかった**——
`autAdjoinRootXPowSubCEquiv_root : σ_η(root) = η • root` が
まさに探していたものだった（次のブロックで訂正）。

**直し方**

`absent` を書く前に、**Mathlib 全体を import した REPL で `#check` を並べる**。

```
mcp__abc3-lean__lean_start  imports: ["Mathlib"]      -- 約 136 秒（1 回だけ）
mcp__abc3-lean__lean_check  #check @Foo / #check @Bar -- 0.02 秒で 10 個
```

★`grep -r` は分単位でかかるうえ、**ディレクトリの見当が外れると空振りする**。
★★`#check` は名前が合っていれば必ず当たり、型まで出る。
★★★`mathlib-gap.json` の `_absentPolicy`（「absent は探索範囲を伴わなければならない」）
の運用として、**探索範囲は「grep したディレクトリ」ではなく
「`#check` した名前の一覧」で書く**のがよい。

## 整数の cast の不等式 —— `exact_mod_cast` が通らないとき

**失敗形**（2026-08-29、`Found/GaloisRep/HtFinJ.lean`）

```lean
-- h : max 0 (-jExp p W) ≤ minDeltaExp p W   (ℤ の不等式)
-- ⊢ ((max 0 (-jExp p W) : ℤ) : ℝ) ≤ ((minDeltaExp p W : ℤ) : ℝ)
exact_mod_cast h
-- mod_cast has type  max 0 (-jExp p W) ≤ minDeltaExp p W
-- but is expected to have type  ↑(max 0 (-jExp p W)) ≤ ↑(minDeltaExp p W)
```

★`max` が挟まると `push_cast` の正規化が両辺で食い違い、`exact_mod_cast` が
「同じものだ」と言えなくなる。

**直し方**

**`Int.cast_le.2 h` を直接使う。** `@[simp, norm_cast] Int.cast_le : (↑m ≤ ↑n) ↔ (m ≤ n)`
なので、cast の向きを手で指定すれば一発で通る。

```lean
exact Int.cast_le.2 (maxJ_le_minDeltaExp p W)
```

★同型の失敗は `Nat.cast_le` / `Rat.cast_le` でも起きる。

## 名前が変わったもの（2026-08-29 に当たったぶん）

| 使えない名前 | 今の名前 |
|---|---|
| `pow_le_pow_left` | **`pow_le_pow_left₀`**（`0 ≤ a → a ≤ b → a^n ≤ b^n`） |
| `div_le_div_iff` | `div_le_div_iff₀` |
| `finsum_le_finsum` | **`finsum_le_finsum'`**（引数は `hf hg h` の順——**台の有限性が先**） |

★`₀` が付く一群は「零因子のない可換モノイド」への一般化で名前が動いたもの。
`exact?` は当たるが遅いので、**まず `₀` を足して試す**。

## `positivity` が `≠ 0` を出せないとき

**失敗形**（2026-08-29、`Found/GenEll/JArchBound.lean`）

```lean
rw [Real.log_mul (by positivity) (by positivity), ...]
-- failed to prove nonzeroness, but it would be possible to prove nonnegativity if desired
```

★`‖jFun τ‖ ≠ 0` は `positivity` には出せない——ノルムは `0` になりうるからである。

**直し方**

**正値性を先に `have` で出しておき、`ne_of_gt` を渡す。**

```lean
have hjpos : (0:ℝ) < ‖jFun τ‖ := by
  rcases (norm_nonneg (jFun τ)).lt_or_eq with h | h
  · exact h
  · rw [hjn, ← h] at hj; simp [Real.log_zero] at hj
rw [Real.log_mul (ne_of_gt hjpos) (by positivity), ...]
```

★`positivity` は「式の形」から出すので、**仮定に隠れている非零性は見えない**。

## `Real.log_two_gt_d9` が「無い」と言われたら —— import が足りないだけ

**失敗形**（2026-08-29、`Found/GaloisRep/Lemma37A.lean`）

Mathlib 全体を import した REPL では通るのに、`lake build` で

```
error: Unknown constant `Real.log_two_gt_d9`
```

★これは**存在しない**のではなく、**そのファイルの import に入っていない**だけである。

**直し方**

`import Mathlib.Analysis.Complex.ExponentialBounds` を足す。
`Real.log_two_gt_d9` / `Real.log_two_lt_d9` / `Real.exp_one_gt_d9` はここにある。

★★**教訓**: REPL（`Mathlib` 全体）で通ったコードをファイルへ移すときは、
**REPL の import とファイルの import の差**を疑う。
「無い」と決めて仮説に逃がす（`hlog : 0.69 ≤ Real.log 2` を受ける）前に、
`grep -rln "<名前>" .lake/packages/mathlib/Mathlib/` で**どのファイルにあるか**を見ること。
★過去に同じ形で仮説へ逃がした記録がある（本ファイルの上のほうの「名前が変わったもの」の節）。

## `WeierstrassCurve.Affine.negY` が Unknown constant になる（2026-08-29、第 593）

REPL（全 Mathlib）では `W.toAffine.negY` が通るのに、ファイルでは

    Invalid field `negY`: The environment does not contain `WeierstrassCurve.negY`
    Unknown constant `WeierstrassCurve.Affine.negY`

になる。`negY` は `Affine/Basic.lean` ではなく **`Affine/Formula.lean`** にある。
`import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Formula` を足す。

★エラーの第 1 行が「`WeierstrassCurve.negY` が無い」と言うのは、
`Affine R := WeierstrassCurve R` が透明なのでドット記法が
`Affine` 名前空間を見つけられず `WeierstrassCurve` に落ちるため。
★★**「不在の定数」ではなく「未 import」を先に疑う**——同じ穴は
`Real.log_two_gt_d9`（`Analysis.Complex.ExponentialBounds`）でも踏んだ。

## 在庫を先に引く——`valAdd_nonneg_iff` を二度書いた（2026-08-29、第 595）

`NeronValuation.lean` に `valAdd_nonneg_iff` を新しく書いたら

    error: `ABC3.Found.GaloisRep.valAdd_nonneg_iff` has already been declared

`DegInf.lean:78` に同名・同 statement が既にあった。
★CLAUDE.md の「在庫」どおり **`node tools/decl-index.mjs` → `.cache/decl-index.txt` を
grep してから書く**。今回は `neronExp` だけ grep して周辺の補題を見ていなかった。
★★下流のファイルに既にある補題は、上流で書き直すと衝突する
（下流で使われている参照は上流の定義に解決されるので、直し方は
「上流の重複を消す」か「新しい定理を下流へ移す」）。

## 三たび——在庫を引かずに書いて名前衝突（2026-08-29、第 612-613）

`Found/GenEll/Uniformization.lean` に `hasDerivAt_weierstrassP` /
`hasDerivAt_derivWeierstrassP` を書いたら

    error: import ABC3.Found.GenEll.WeierstrassODE failed, environment already
    contains 'ABC3.Found.GenEll.hasDerivAt_weierstrassP' from ...Uniformization

★★**同じ名前空間の別ファイルに既にあった**（`Found/GenEll/WeierstrassODE.lean:62`）。
しかも同ファイルには **`deriv_derivWeierstrassP : deriv ℘′ z = 6·℘(z)² − g₂/2`**
（`℘` の 2 階微分）まであり、これは今まさに必要としていたものだった。

★CLAUDE.md の「在庫」を守ること——**書く前に `node tools/decl-index.mjs` →
`.cache/decl-index.txt` を grep する**。名前だけでなく**周辺の定理**も見る。
☆本セッションで 2 度目（第 595 の `valAdd_nonneg_iff`）。
★★★このエラーは `lake build` の import 段で出るので、
**単一ファイルのビルドでは気づけない**——`lake build ABC3` を必ず通すこと。


## `set` で入れた局所定義を `ring` が展開してしまう（2026-08-29、第 666）

`set ω₁' := η₁ / (l : ℂ) with hω₁'` としたあと `linear_combination` を書いたら

    ring failed, ring expressions not equal
    ⊢ P.ω₁ * 2 - P.ω₁ * ↑p * ↑a₁ * ↑l * (↑l)⁻¹ * 2 - … = 0

★**`set` は `let` 束縛を作るので `ring` / `linear_combination` が中身まで
展開する**（`(↑l)⁻¹` が出てきたのがその証拠）。抽象的な等式
（`(l : ℂ) * ω₁' = η₁`）を先に `have` で取ってから

    clear_value ω₁'

で本体を落とすと、以降 `ω₁'` は不透明な局所定数として扱われ、
`linear_combination (係数) * hlω'` が期待通り効く。

☆同じ理由で **`rintro w (rfl | rfl)` は `set` 変数を消してしまう**ことがある
（`w = ω₁'` の `rfl` が `ω₁'` の側を除去して `Unknown identifier ω₁'`）。
`intro w hw` → `simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hw`
→ `rcases hw with h1 | h1` → `rw [h1]` と書けば向きを固定できる。

## `Nat.mul_le_mul_right` の引数順（2026-08-29、第 667）

    calc n = 1 * n := (Nat.one_mul n).symm
      _ ≤ g * n := Nat.mul_le_mul_right n hgpos   -- ← 期待と違う辺が出る
      _ = l := hn.symm                            -- Type mismatch: g * n = n

★`Nat.mul_le_mul_right` は版によって明示引数の位置が違う。
不等式は **`Nat.le_mul_of_pos_left n hgpos : n ≤ g * n`** を使うほうが安全。
☆`l = g * l` から `g = 1` を出すのも `rw [← hnl]` だと `l` を全部書き換えて
しまうので、`conv_lhs => rw [hn]` で片側だけ触ってから `rw [hnl]` とする。


## `HasDerivAt` の合成補題は Pi 形の関数を作る（2026-08-29、第 673）

    have hb := ((h1.pow 2).const_mul (6 : ℂ)).sub_const (P.g₂ / 2)
    rw [hb.deriv]   -- ← 失敗

    Tactic `rewrite` failed: Did not find an occurrence of the pattern
      deriv (fun x ↦ 6 * (℘[P] ^ 2) x - P.g₂ / 2) w
    in the target expression
      deriv (fun z ↦ 6 * ℘[P] z ^ 2 - P.g₂ / 2) w

★`HasDerivAt.pow` は `℘ ^ 2`（`Pi.pow`）を、`HasDerivAt.mul` は `℘ * ℘'`
（`Pi.mul`）を作る。ラムダ形の目標とは**構文が違う**。直し方:

    have h2 : HasDerivAt (fun z : ℂ => 6 * P.weierstrassP z ^ 2 - P.g₂ / 2) _ w :=
      hb.congr_of_eventuallyEq (by filter_upwards with z; simp only [Pi.pow_apply])

☆導関数の値は `_` にしておけば `hb` から推論される。
☆`Finset.analyticAt_sum` も同じ罠で、結論は `AnalyticAt 𝕜 (∑ n ∈ N, f n) c`
（Pi 和）なので `fun z => ∑ ...` には `.congr` ＋ `simp [Finset.sum_apply]` が要る。

## 対称な目標で `rw [h1, h2]` が両辺を書き換える（2026-08-29、第 670）

    have hSS : S = -S := by rw [h1, h2]
    -- ⊢ -∑ ... = - -∑ ...

★目標 `∑ = -∑` の**両辺に**同じ部分項があると `rw` は両方を書き換える。
`h1.trans h2` と項で書くか、`conv_lhs => rw [...]` で片側に閉じ込める。


## `Point` の群構造は `[DecidableEq F]` を要求する（2026-08-29、第 686）

    failed to synthesize instance of type class
      AddMonoid W.toAffine.Point

★mathlib の `WeierstrassCurve.Affine.Point.instAddCommGroup` は
`variable [DecidableEq F]` の節にある。`+` だけなら `instance : Add W.Point`
（`DecidableEq` 不要）で足りるが、**`n • ` や `addOrderOf` を使うと `AddMonoid` が要り、
そこで詰まる**。`open scoped Classical in` を定理の直前に置けば通る。
☆`[W.IsElliptic]` を足しても直らない——足りないのは `DecidableEq` のほうである。

## `u⁻¹` の掛け算は `ring` では消えない（2026-08-29、第 689）

`Units` の `u` と `u⁻¹` は `ring` にとって別々の不定元なので、
`u⁻²·u² · x = x` は `ring` で閉じない。`hu : (u:F) * (u⁻¹:F) = 1`（`C.u.mul_inv`）を
**係数つきで** `linear_combination` に渡す:

    linear_combination (x' * ((C.u : F) * (C.u⁻¹ : F) + 1)) * hu

★`(uu⁻¹)² − 1 = (uu⁻¹ − 1)(uu⁻¹ + 1)` の因数分解が係数の形を決める。
3 乗なら `(uu⁻¹)² + uu⁻¹ + 1` が係数になる。
☆逆に `variableChange_a*` は `u⁻¹` **だけ**で書かれているので、そこは素の `ring` で閉じる。


## 在庫: `IsDedekindDomain.HeightOneSpectrum.under` は mathlib にある

**失敗形**: `P ∩ 𝓞L` を `HeightOneSpectrum` として作る定義を自分で書いた。

**実際**: `Mathlib/RingTheory/DedekindDomain/Ideal/Lemmas.lean` に
`IsDedekindDomain.HeightOneSpectrum.under (A) (w : HeightOneSpectrum B) : HeightOneSpectrum A`
がある(`@[simps]` 付き)。`HeightOneSpectrum.comap` が全射環準同型限定なのに気を取られて
`under` を見落とした。

**直し方**: `decl-index.txt` を `HeightOneSpectrum` で grep するとき、`comap` だけでなく
`under`・`over`・`liesOver` も見る。イデアルの contraction は mathlib では `Ideal.under`
という名前である(`comap` ではない)。

## 在庫: `finite_j_of_htFalt_le` は既にある（Northcott は済んでいた）

**失敗形**: 「`ht^Falt` が有界な類は有限」（Northcott）が無いと思い込み、
`Skeleton/GenEll/NorthcottJ.lean` に節点を立てた（第 743）。

**実際**: `Found/GaloisRep/NorthcottHtJ.lean` の `finite_j_of_htFalt_le`（`§9-1005`）が
**無条件で証明済み**だった。`Found/GenEll/NorthcottImage.lean`（`§9-950`）が
`htJ` の Northcott 性を与えている。★mathlib に instance が無いことは正しかったが、
**本プロジェクトが自前で作っていた**。

**直し方**: 新しい節点を立てる前に、まず `node tools/decl-index.mjs` で
`.cache/decl-index.txt` を作り、`finite`・`Northcott`・`Finite` で grep する。
★★特に「§n-xxxx で済んでいる」と他のファイルの docstring が書いていないか、
`grep -rn "無条件" ABC3/Found/` で見る。

**副産物**: この測定で `SSCurve.fld : Type` ＋ 埋め込みという設計が
`finite_j_of_htFalt_le`（族の定義体を `IntermediateField ℚ ℂ` で受ける）と
噛み合わないことも分かり、`SSCurve.K : IntermediateField ℚ ℂ` に直した（第 753）。

## 数体の判別式・円分体（2026-08-31、第 780-784 で使った mathlib の在庫）

「`l` が `L` で不分岐」を **`¬ (l : ℤ) ∣ NumberField.discr L`** として持つと、
mathlib の以下がそのまま噛み合う。分岐理論の API を探す必要はない。

* `NumberField.discr_dvd_discr : discr K ∣ discr L`（`K ⊆ L`）
* `NumberField.not_dvd_discr_iff_forall_mem : ¬ p ∣ discr K ↔ ∀ P ∋ p, IsUnramifiedAt`
* `NumberField.abs_discr_gt_two : 1 < finrank ℚ K → 2 < |discr K|`
* `NumberField.finrank_eq_one_of_unramified`（至る所不分岐な数体は `ℚ`）
* `NumberField.linearDisjoint_of_isGalois_isCoprime_discr`
  —— **`K₁/ℚ` が Galois で `disc` が互いに素なら線型無関連**（これが要）
* `IsCyclotomicExtension.Rat.discr_prime_pow : NumberField.discr K = ±p^m`
* `IsCyclotomicExtension.Rat.finrank : finrank ℚ ℚ(ζ_n) = n.totient`
* `IntermediateField.LinearDisjoint.adjoin_rank_eq_rank_left_of_isAlgebraic`
  —— `[L(A):L] = [A:F]`

### 失敗形と直し方

* `Int.coe_nat_prime` / `Int.natCast_prime` は**無い**。
  `Prime (l : ℤ)` は `rw [Int.prime_iff_natAbs_prime]; simpa using (Fact.out : l.Prime)`。
* `IsCyclotomicExtension.isGalois` の第 1 引数は **`Set ℕ`**。`{l ^ k}` と書く（`(l ^ k)` は型エラー）。
* `IsPrimitiveRoot.powerBasis K` は **`[NeZero ((n : ℕ) : K)]`** を要る。
  `⟨(Nat.cast_ne_zero (R := K)).2 (NeZero.ne n)⟩` で供給する。
* `IsPrimitiveRoot.minpoly_eq_cyclotomic_of_irreducible` の向きは
  **`cyclotomic n K = minpoly K μ`**（左右が直感と逆）。`.symm.trans` で使う。
* `IsPrimitiveRoot.minpoly_dvd_cyclotomic` は `μ` が **`K` 自身**にある場合のみ。
  拡大体の `μ` には `minpoly.dvd K μ (aeval μ (cyclotomic n K) = 0)` を使い、
  根であることは `rw [aeval_def, ← eval_map, map_cyclotomic]` で `isRoot_cyclotomic` に落とす。
* `Normal.of_isAlgClosed` は**無い**。`IsAlgClosure.normal B Ω`（`IsAlgClosure` を
  `⟨inferInstance, inferInstance⟩` で作ってから）。
* 既約性を環同型で移すのは `MulEquiv.irreducible_iff (Polynomial.mapEquiv e)`
  （`f` は**明示引数**なので `.irreducible_iff` のドット記法は効かない）。
* `IntermediateField.restrict_algEquiv` の向きは `↥E ≃ₐ[F] ↥(E.restrict h)`。
* `le_sup_right` を包含として使うときは型注釈が要る: `(le_sup_right : Z ≤ M) hζZ`。
* `NumberField` のインスタンスは `⟨⟩`（明示欄が無い）。`Module.Finite.trans (R := ℚ) (A := B) (M := M)`
  は引数を取らない（`inferInstance` を渡すと "Function expected"）。

## 配管：新規ファイルの登録先（2026-08-31、第 800）

`lean/ABC3.lean` は **アグリゲータだけ**を import する（`ABC3.Meta.Claim` /
`ABC3.Interface` / `ABC3.Skeleton` / `ABC3.Found` / `ABC3.Gap` / `ABC3.Check`）。
★新しい `.lean` を作ったら **`lean/ABC3/Found.lean`（等）に `import` を足す**こと。
`ABC3.lean` に足しても行が一致せず sed が黙って何もしない。

☆症状: 個別ビルド（`lake build ABC3.Found.GaloisRep.Foo`）は通り、
`node tools/check.mjs` も PASS するのに、`lake build ABC3` がそのファイルを
一度もコンパイルしていない（他のファイルから import されていなければ）。

## Tate の `I` 進級数：adic 添字と `q` 次数はずれる（2026-08-31、第 818）

`tateXtail (w, q) = adicSum (n ↦ q^n · ∑_{d∣n} d·w^d)` である。
★`w = ζ` のときは `q` 次数 = adic 添字 `n` だが、
**`w = q·ζ^{-1}` のときは `w^d = q^d ζ^{-d}` なので `q` 次数は `n + d`** になる。

☆したがって「adic 添字ごとに係数を比べる」形の照合は、
先に **`q` 次数に揃えた形**（古典形 `X(u) = u/(1−u)² + ∑_N (∑_{d∣N} d(u^d + u^{−d} − 2)) q^N`）
へ直してからでないとできない。

★道具は `AdicFubini.lean` の `adicSum_reindex_mul`・`adicSum_fubini`、
`AdicMul.lean` の `adicSum_mul`。

## IsLocalization.lift が通らない——インスタンスの菱形（2026-08-31、第 845）

```
univDual has type @RingHom TateBase (TrivSqZeroExt ..) AddMonoidAlgebra.nonAssocSemiring
                                                       nonAssocSemiring
but is expected  @RingHom TateBase ?P  AddMonoidAlgebra.commSemiring.toNonAssocSemiring
                                                       CommSemiring.toSemiring.toNonAssocSemiring
```

`TrivSqZeroExt` も `AddMonoidAlgebra` も `nonAssocSemiring` を**独立に**宣言しており、
`CommSemiring` 経由の道と構文的に一致しない。`by exact` も `(g := ..)` も効かない。

**直し方（両方使う）**

* **標的側**——`def MyType := TrivSqZeroExt R M`（`abbrev` ではなく `def`）で包み、
  `noncomputable instance : CommRing MyType := inferInstanceAs (CommRing (TrivSqZeroExt R M))`
  を 1 つだけ与える。`def` は reducible でないので探索は唯一の道を通る。
* **始域側**——準同型を `MvPolynomial.eval₂Hom` で作る。
  戻り値の型が `[CommSemiring R]` 由来になるので `IsLocalization.lift` と揃う。

包んだ型の上で `TrivSqZeroExt.snd_mul` などを使うには、まず
`theorem mul_eq (x y) : (show TrivSqZeroExt R M from x * y) = (show .. from x) * (show .. from y) := rfl`
を置き、`simp only [eps, re, mul_eq, TrivSqZeroExt.snd_mul]` の順で展開する。
☆`simpa … using h` は `h` を `True` にしてしまうことがあるので `simp only` で順に展開する。

☆付随して見えた失敗形:

* `Int.induction_on` の case 名は `zero` / `succ` / `pred`（`hz`/`hp`/`hn` ではない）
* `MvPolynomial.induction_on` の case 名は `C` / `add` / `mul_X`（`h_C` 等ではない）
* `Fin 2` の `fin_cases i` は `X ((fun i ↦ i) ⟨0, ⋯⟩)` を作り `rw [… X 0]` が失敗する。
  `theorem foo : ∀ i : Fin 2, … | 0 => … | 1 => …` の**等式コンパイラ**で書く。
* `IsLocalization.induction_on` は無い。`IsLocalization.surj` を使い
  `obtain ⟨⟨a, s⟩, hz⟩ := IsLocalization.surj (M := …) x`（`hz : x * ι s = ι a`）とする。
* 局所化への**微分の延長**は mathlib に無い（`Mathlib/RingTheory/Derivation/` を
  `Localization` で grep して 0 件）。双対数 + `IsLocalization.lift` で自分で作る。

## Python で Lean ファイルを書き換えるときは一時ファイル経由で（2026-08-31、第 872）

```python
io.open(p, "w", encoding="utf-8").write(s)   # ✕ 危険
```

`io.open(..., "w")` は**開いた瞬間に切り詰める**ので、`write` が例外で落ちると
**原本が空になる**。実際にサロゲート対を書こうとして
`UnicodeEncodeError: surrogates not allowed` で落ち、204 行のファイルが 0 行になった。

```python
tmp = p + ".tmp"
with io.open(tmp, "w", encoding="utf-8", newline="\n") as f:
    f.write(out)
os.replace(tmp, p)          # ○ 途中で落ちても原本は無事
```

☆アストラル面の文字は **8 桁の** `\U0001D4DE` のように書く。
6 桁の `\ud835` + `\udcde` のようなサロゲート対は Python 3 では書けない。
☆似た形の字に注意——本プロジェクトの整数環は `𝓞`（\U0001D4DE）であり、`𝒪`（\U0001D4AA）ではない。
☆大きな書き換えの前に commit しておけば、万一のとき `git checkout -- <path>` で戻せる。
☆在庫の確認を先に——`jExp_congr_j` は既に `Found/GaloisRep/HtFaltJ.lean` にあった。

## Python の「次の `def` まで消す」スプライスが定理を巻き込む（2026-08-31、第 873→881）

不要になった `def foo.needs_old` を消そうとして、こう書いた:

```python
i = s.index("def tateModel_of_quot_mu.needs_old")
j = s.index("\ndef ", i + 10)          # ★次の def まで
out = s[:i] + s[j+1:]
```

`needs_old` ブロックの直後に `def` が来ない（`theorem` が 3 つ挟まる）と、
`s.index("\ndef ", ...)` は**それらを飛び越えた先**を指す。結果
`theorem c4_velu_tate`・`c6_velu_tate`・`j_velu_tate_mu` が丸ごと消えた。

★**気付けなかった理由**——消えた 3 つを名前で参照している所が無く
（`.needs` の中の**文字列**でしか引かれていない）、`lake build` も
`check.mjs` も通ってしまった。8 ブロック後に `j_velu_tate_mu_map` を
書いて初めて `Unknown identifier 'c4_velu_tate'` で露見した。

☆直し方は 2 つ:

1. 終端を**次の `def` ではなく `\n\n`**（空行）にする
2. もっと良いのは、削る範囲を**両端の文字列で挟んで指定**する:

```python
i = s.index("def foo.needs_old")
j = s.index("def bar.src", i)          # ★次に残したいものを名指しする
out = s[:i] + s[j:]
```

★そして**削った後に必ず** `grep -c "^theorem " file` の前後を比べる。
行数が減っているのは当たり前なので、**宣言の数**を見ること。

## `IsElliptic` のインスタンスが `rw` の motive を壊す（2026-08-31、第 882・890）

`W.j` は `[W.IsElliptic]` を要求する。したがって `h : A = B`（曲線の等式）で
`rw [h]` すると motive が `fun x => x.j = …` になり、`x` に対して
`IsElliptic` が合成できず **motive is not type correct** で落ちる。
同じことが `tatePtPair a w q hq (haw : a*w=q) …`（証明を引数に取る）でも起きる。

☆逃げ道は 3 つある:

1. **`j` を経由しない**。`c4_cube_mul_Delta_of_j_eq` で分母を払い、
   `Δ ↦ u⁻¹¹²Δ`・`c₄ ↦ u⁻¹⁴c₄` を使って**積の形**で計算する
   （`Found/GaloisRep/TateParamJ.lean`）。変数変換の `u` は両辺で同じ重みを拾って消える。
2. **`subst` で潰す**。等式の片側が**局所変数**なら `subst` は motive を作らないので通る。
   `tateParam_quot_velu`（第 891）は `hW' : W' = …` を `subst` してから `rfl` で済ませている。
3. **`_congr` 補題を作る**。`tatePtPair_congr`（第 884）は値の引数だけを
   `subst` で潰し、証明の引数は暗黙にしてある。

## `open scoped Classical` と `[DecidableEq K]` を同時に置くと `Finset.image` が食い違う（2026-08-31、第 890）

`Finset.image` は `DecidableEq` をデータとして持つ。ファイル A で
`open scoped Classical` だけ、ファイル B で `[DecidableEq K]` も宣言していると、
B の目標の `image` は `inst✝` を、A の補題は `Classical.propDecidable` を持ち、
`exact` が **Type mismatch** で落ちる（型は「同じに見える」ので原因が分かりにくい）。

☆直し方は**どちらかに揃える**こと。本プロジェクトは `open scoped Classical` 側に
揃えた（`[DecidableEq K]` を variable から外す）。
★`tatePhi` のように `[DecidableEq K]` を要求する定義も、Classical があれば通る。

## `have i : C X := inferInstance` を並べるとインスタンスが壊れる（2026-08-31、第 898）

「必要なインスタンスが全部あるか」を 1 つの `example` の中で

```lean
example ... : True := by
  have i1 : Field Lv := inferInstance
  have i2 : Algebra L Lv := inferInstance   -- ★ここで落ちる
  ...
```

と並べて確かめてはいけない。`have i1` を置いた瞬間、局所文脈に
**2 つ目の `Field Lv`** が入る（大域のインスタンスと `i1`）。
`Algebra L Lv` は `Semiring Lv` を経由するので、どちらの `Field` から
降りるかで菱形になり **failed to synthesize** で落ちる。

★紛らわしいのは、同じ import・同じ `open` でも
**`i2` だけを書いたファイルは通る**ことである（原因が import に見える）。

☆直し方: **各インスタンスを別々の `example` で確かめる**。

```lean
example ... : Algebra L (p.adicCompletion L) := inferInstance
example ... : IsFractionRing (p.adicCompletionIntegers L) (p.adicCompletion L) := inferInstance
```

★`Algebra` のようにデータを持つクラスは `noncomputable example` にすること
（`UniformSpace.Completion.instField` 等が noncomputable なので）。

## `omit` / `open ... in` は docstring の**前**に置く（2026-09-01、第 932・940）

宣言修飾子（`omit [C] in`、`open scoped Classical in`、`set_option ... in`）を
docstring と宣言のあいだに挟むと

```
error: unexpected token 'omit'; expected 'lemma'
```

になる。docstring は宣言の**直前**でなければならない。

```lean
-- ★正しい
omit [CharZero K] in
/-- ★説明 -/
theorem foo : ... := ...

-- ✗ 落ちる
/-- ★説明 -/
omit [CharZero K] in
theorem foo : ... := ...
```

☆`section variable` に入れた不要なインスタンスは、この `omit` で個別に外す。
★`variable` を分けるより、`omit` 1 行のほうが差分が小さい。

## `′`(U+2032 PRIME) は Lean の識別子に使えない

`variable {R′ : Type}` は `expected token` で落ちる。docstring や文字列の中では
問題ないので、**識別子だけ ASCII の `'` にする**（`R'`・`h'`・`W'`）。
第 944 で 79 箇所を一括置換した。

## `IsElliptic` の motive 壊れ、三度目——`X.j` の `X` を書き換えるとき

`rw [hEq]`（`hEq : X = Y`）を `... = X.j` の上でやると
`motive is not type correct`（`j` が `X.IsElliptic` を暗黙に取るため）。
★`ABC3.Found.GenEll.j_congr_curve hEq : X.j = Y.j` を **simp 補題として渡す**。
`rw` だと後続の `map_j` の出現位置がずれるので `simp only [...]` が安全。

## 「全称で受けたデータ」は充足不能になりうる

第 948 で `hv : ∀ ζ, IsPrimitiveRoot ζ l → v = ∑ …(ζ)` と書いたが、これは
**1 つの `v` が全ての原始根について成り立つ**ことを要求しており、和が `ζ`
に依らないことを別に証明しない限り満たせない。★正しくは

    `hvw : ∀ ζ, IsPrimitiveRoot ζ l → ∃ v w, … ∧ … ∧ …`

と**存在量化を内側に入れる**（第 952 で直した）。
☆補助データ（ここでは `v`・`w`）が結論に現れないなら、必ず内側の ∃ にする。

## 在庫を引く前に書き始めない（第 958 で再発）

`veluU_negY`・`veluGy_negY` を新規に書こうとしたら
`has already been declared` で落ちた——`Found/GenEll/Velu.lean:740` に既にあった。
★CLAUDE.md の「在庫」の通り、**書く前に `node tools/decl-index.mjs` →
`.cache/decl-index.txt` を grep する**。名前が思いつく補題ほど既にある。

## `Point` の `•` は `open scoped Classical` の側で取れている（第 966）

`rhPoint_nsmul` などの `n • Pt` は、宣言側のファイルが `open scoped Classical`
なので **Classical の `DecidableEq` から来る加法群**の `SMul` を使っている。
★呼ぶ側で `[DecidableEq F]` を束縛すると別インスタンスになり、`rw` が
「パターンが見つからない」で落ちる。
☆呼ぶ側も `open scoped Classical in` にして `DecidableEq` は束縛しないこと。

## 在庫の引き忘れ、二度目（第 967）

`addOrderOf_rhPoint`・`addOrderOf_vcPoint` を新規に書こうとして
`has already been declared`。**`Found/GenEll/PointVariableChange.lean:452, 1001`**
に既にあった。★1 セッションで 2 回同じ穴に落ちた（第 958 と第 967）。
☆対策: 補題を書く前に、**その回に使う名前を全部まとめて** decl-index に
grep する（1 つずつ思い出した順に引くと漏れる）。

## 同名の補題が 2 つの名前空間にある（第 970）

`vcPoint` は `ABC3.Found.GaloisRep`（`W C` の順）と `ABC3.Found.GenEll`（`C W` の順）の
**両方にある**。★`namespace ABC3.Found.GaloisRep` の中で裸の `vcPoint` と書くと
GaloisRep 側が選ばれ、`Application type mismatch` になる。
☆MCP の `lean_check` は top-level ＋ `open` で試すので**この衝突が出ない**。
ファイルに移すときは名前空間を明示すること。

## `Valued.mem_nhds` は `restrict` と `ValueGroup₀` の言葉（第 989）

`s ∈ 𝓝 x` を**作る**とき、`rw [Valued.mem_nhds]` の後の目標は
`Valued.v.restrict (z − x) < ↑γ` の形になる。
★`γ := 1` を渡し、`Valuation.restrict_lt_iff_lt_embedding` で書き換えてから
`simp only [Units.val_one, map_one]` で `Valued.v (z − x) < 1` に落とす。
☆第 897・943 でも同じ場所で止まっている——3 度目。

## `rw` の後に `X = X` が残る(第 999)

**失敗形**: `rw [hPz]` で両辺が同じ形になったのに「unsolved goals」で
`⊢ X = X`(表示上は完全に同一)が残る。

**なぜ**: `rw` の自動 `rfl` は reducible 透明度で走る。
`Finset.image` の `DecidableEq` インスタンスや `tateCurveAt` に渡した証明項が
片側だけ別経路で作られていると、表示は同じでも reducible には合わない。

**直し方**: `rw [hPz]` の直後に **`rfl` を明示的に書く**。
`rfl` タクティクは default 透明度＋証明無関係性を使うので通る。

```lean
have hcurveEq : veluQuotientFull W (image (fun k => pointCoords (k • P)) s)
    = veluQuotientFull W (image (fun k => pointCoords (k • tatePhi S hΔ c)) s) := by
  rw [hPz]
  rfl        -- ★これが要る
```

**併せて**: `j` を跨ぐ書き換えは motive が壊れる(`IsElliptic` が邪魔)。
**曲線の水準で等式を作ってから** `ABC3.Found.GenEll.j_congr_curve` で `j` に移すと安全。

## Python でファイルを書き換えるとき、開いた瞬間に空になる(第 1013)

**失敗形**: `io.open(path,'w')` に渡す文字列の組み立てで例外が出ると、
**ファイルはすでに truncate されていて 0 バイトになる**。
(実例: `u'𝔪'` のようなサロゲート対を書くと
`UnicodeEncodeError: surrogates not allowed` が出る。𝔪 は `\U0001D52A`。)

**直し方**: **書き切ってから置換する**。

```python
tmp = src + '.tmp'
io.open(tmp, 'w', encoding='utf-8', newline='\n').write(out)
os.replace(tmp, src)
```

**併せて**: Lean のコード片は Write ツールで別ファイルに書き、
Python はそれを読んで挿入するだけにすると、エスケープ事故そのものが起きない。
壊したときは `git checkout HEAD -- <path>` で戻す(だからこまめに commit する)。

## monic・次数の証明は `monicity!` / `compute_degree!`(第 1031)

`X^2 + C a * X - C b` のような多項式について
`Monic` と `natDegree = 2` を手で示すと `degree_sub_le` / `max_le_iff` /
`degree_C_mul_le` を並べることになり、`WithBot ℕ` と `ℕ` の往復で崩れやすい。

**mathlib のタクティクを使う**:

```lean
theorem monic_p : p.Monic := by rw [p]; monicity!
theorem natDegree_p : p.natDegree = 2 := by rw [p]; compute_degree!
```

`!` 付きは残った副目標を `norm_num`/`assumption` で閉じにいく。

## `Ideal.Quotient.mk (maximalIdeal R)` は `Field` 探索が爆発する(第 1036)

**失敗形**: `Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R)))` を
体上の補題(`[Field k]` を要求)に当てると
`(deterministic) timeout at whnf` になる。
`Field (R ⧸ 𝔪)` を見つけるのに `ResidueField R` まで展開する必要があるため。

**直し方**: `show` で **`IsLocalRing.residue R`** の形に言い直してから当てる。
`residue R : R →+* IsLocalRing.ResidueField R` なので `Field` が即座に見つかる。

```lean
  show Irreducible (f.map (IsLocalRing.residue R))
  refine some_lemma (hf.map (IsLocalRing.residue R)) ?_ hns   -- 6.9 秒 → 0.04 秒
```

**併せて**: `refine g h ?_ ?_` で第 1 引数の暗黙変数がゴールから決まらないときは
`(q := ...)` で明示する。それでも遅いときは、**大きな項を含む補題を
小さい文脈(`Found/` の抽象的な補題)に切り出す**と桁で速くなる。

## `↑` と `algebraMap` が一致せず `rw` が刺さらない(第 1040)

**失敗形**: `HeightOneSpectrum.valuation_of_algebraMap` を `rw [←]` すると
ゴールに `↑↑m` が現れるが、これは `algebraMap ... ↑m` と**構文的に別物**で、
`map_natCast` も自作の `have e1 : algebraMap ... = ...` も刺さらない。
(`𝓞 L → L` の coe と `algebraMap (𝓞 L) L` の食い違い。)

**直し方**: 等式を `rw` で押し込まず、**`congr 1` に落として `push_cast` で潰す**。

```lean
  rw [← HeightOneSpectrum.valuation_of_algebraMap (K := Lv), 
      ← HeightOneSpectrum.valuation_of_algebraMap (K := L)]
  rw [← hp (algebraMap (𝓞 L) L ((m : ℕ) : 𝓞 L))]   -- hp の側を先に合わせる
  congr 1
  push_cast
  ring
```

**併せて**: `set_option linter.tacticCheckInstances true` を勧めるノートが出たら、
インスタンスのダイヤモンド由来なので `rw` ではなく `congr`／`convert` に切り替える。

## `#no-guard` は行末までのコメント —— 同じ行の後続コマンドが消える

ABC3 の Bash ガードを外す `#no-guard` は**シェルのコメント**なので、`#` から
行末までが丸ごと捨てられる。次は前半しか実行されない:

```bash
ls some/dir #no-guard; echo "=== 次 ==="; grep -rn foo some/dir   # ← echo も grep も実行されない
```

出力に「後半の echo が出ていない」ことでしか気づけないため、
**測定結果を「0 件だった」と誤読する**（第 1068 → 第 1069 で訂正）。

直し方: `#no-guard` は**コマンド全体の最後**、できれば独立した最終行に置く。

```bash
ls some/dir
echo "=== 次 ==="
grep -rn foo some/dir
#no-guard
```

☆教訓: grep が「0 件」を返したときは、**同じコマンド内の直前の echo が
出力されているか**を必ず確かめる。出ていなければ grep は走っていない。

## 在庫は「mathlib に無い」で止めず、**プロジェクト内も**引く

第 1067-1069 で「楕円曲線の形式群は mathlib に無い」「分点多項式と捻れ点の橋が無い」と
測定したが、**プロジェクト内には既に**

- `redPoint` / `redHom` / `redPoint_add` / `redPoint_nsmul`（`Found/GaloisRep/RedKernel.lean`）
- `one_lt_val_addX_of_infinity`（`E₁` が加法で閉じている、`InfinityKer.lean`）
- `val_y_sq_eq_val_x_cube`（`w(y)² = w(x)³`、`Infinity.lean`）
- `Psi` / `Phi` / `PSq` / `PsiRec` / `PsiDouble` / `OmegaMulPsi` / `Eds*`

があった（Theorem 3.8 の Weil 対のために積まれていた）。

☆教訓: `node tools/decl-index.mjs` の grep は**mathlib を調べる前に**やる。
★特に「別の定理のために積んだ機械」は名前が違うので、
概念名（`formalGroup`）ではなく**振る舞い**（`red`, `val_.*x`, `_add`）で引く。

## Bash のヒアドキュメントは**バックスラッシュを 1 段潰す**

`<<'EOF'`（クォートつき）でも、このセッションの Bash 経路では
`\` が `\` に、`/\/g` が `/\/g` に落ちる。JavaScript の正規表現や
Python の `'\'` を含むスクリプトをヒアドキュメントで書くと**壊れる**。

```
# 壊れる
console.log(p.replace(/\/g, '/'))   → /\/g  で SyntaxError

# 壊れない(バックスラッシュを一切書かない)
console.log(p.split(path.sep).join('/'))
```

☆教訓: バックスラッシュを含むファイル内容は **Write / Edit ツールで書く**。
★ヒアドキュメントは「バックスラッシュが 1 個も無い」ときだけ使う。

## 日本語の `.src` / docstring が**文字化けする**ことがある

第 1119 前後で `Lemma 3.5(Vélu の v の…)` が
`Lemma 3.5(VÃ©lu ã® v ã®…)` になっていた（`MuPairDenomFree.lean` 1 行、
`MuDenomFreeSum.lean` 6 行、`AdicEvalGen.lean` 1 行）。
**UTF-8 のバイト列を latin1 として読んで書き戻した**形である。
`lake build` は通る（コメント・文字列リテラルなので型に影響しない）ので
ビルドでは捕まらない。

```bash
node tools/mojibake.mjs        # 走査(検出したら終了コード 1)
node tools/mojibake.mjs --fix  # 復元して書き戻す
```

機構: 化けた行は「latin1 で符号化 → UTF-8 で復号」が成功して**別の文字列**になる。
正常な日本語行は U+00FF を超える文字を含むので latin1 で符号化できず素通りする。

☆教訓: 日本語を含む行をスクリプトで書いたら `node tools/mojibake.mjs` を回す。
★`.src` の `item` が化けると進捗指標が原典の項目名と照合できなくなる。

## Python で Lean を書き換えるときの `𝓞`（サロゲート）

**失敗形**: heredoc の Python に `'\ud835\udcde'` と書くと、Python 3 では
孤立サロゲート 2 文字になり、ファイル中の `𝓞` と一致しない。
`assert old in s` が落ちる（か、書き込むと壊れたファイルになる）。

**直し方**: BMP 外の文字は `'\U0001D4DE'`（大文字 `U` + 8 桁）で書く。
`𝓞` = `\U0001D4DE`、`𝓞 L` は `HeightOneSpectrum (𝓞 L)` に現れる。

**もう一つ**: heredoc の中に長い Python を入れると bash が
`unexpected EOF while looking for matching` で落ちることがある。
そのときは Write ツールで `.py` を作ってから実行する。

## `Units.ext` のあとは `push_cast` ではなく `simp only [Units.val_mul, Units.val_mk0]`

**失敗形**: `Units.mk0 (a^2*x) h = Units.mk0 a ha * Units.mk0 a ha * Units.mk0 x hx`
を `apply Units.ext; push_cast; ring` で閉じようとすると、`push_cast` が
右辺だけを `↑a^2 * ↑x` にまとめて左辺の `↑(Units.mk0 …)` を残し、`ring` が閉じない。

**直し方**: `apply Units.ext; simp only [Units.val_mul, Units.val_mk0]; ring`。

## `Affine.Point` の `DecidableEq` は `open scoped Classical` では揃わない

**失敗形**: `WeierstrassCurve.Affine.Point` の `AddCommGroup` は `[DecidableEq F]` を取る。
一般の体 `L` を量化した命題（`Lemma35Unconditional` 等）は
`open scoped Classical` の下で `fun a b => Classical.propDecidable (a = b)` を拾っているが、
具体型（`↑E.K`、すなわち `ℂ` の中間体）では `Subtype.instDecidableEq` の方が優先される。
結果、`addOrderOf Q = l` が**別の群構造の上の命題**になり、

    Application type mismatch: ... @Affine.Point.instAddCommGroup ... fun a b ↦ a.instDecidableEq b
                               ... 期待されるのは ... fun a b ↦ Classical.propDecidable (a = b)

が出る。☆そのあと `isDefEq` / `whnf` のタイムアウトが続くこともある。

**直し方**: 定義と証明の両方で `letI` で固定する。

```lean
letI : DecidableEq E.fld := fun a b => Classical.propDecidable (a = b)
```

★`open scoped Classical` だけでは**具体型に対しては効かない**。
☆`set E' := … with h` が仮説を書き換えてくれないのも同じ原因——
項が構文的に一致していない。`set` をやめて式を直接渡し、`hE' := rfl` で済ませる。

## `Finset.sum_image` は高階単一化に失敗する

**失敗形**: `(Finset.sum_image hinj).symm` を `calc` の一段に直接置くと

    Type mismatch: ∑ x ∈ S, ?m (semiPair Φ x) = ∑ x ∈ Finset.image ..., ?m x
    but is expected to have type ∑ z ∈ S, veluV2 W (semiPair Φ z).1 ...

——`?m`（和を取る関数）が推論できない。

**直し方**: 型を明示した `have` に置いてから使う。

```lean
have hsum : (∑ z ∈ S.image f, g z) = ∑ z ∈ S, g (f z) := Finset.sum_image hinj
rw [himg] at hsum
```

**もう一つ**: 仮説に `S.image f = S` を置くと **statement 側で `DecidableEq` が要る**
（`classical` は証明の中だけ）。★仮説は `∀ z ∈ S, f z ∈ S`（安定性）にして、
像の等式は `Finset.eq_of_subset_of_card_le` で証明内に作る。

## `open scoped Classical in` も docstring の**前**に置く

**失敗形**: docstring のあとに `open scoped Classical in` を置くと

    unexpected token 'open'; expected 'lemma'

★`set_option … in` とまったく同じ穴である。

**直し方**: `open scoped Classical in` → `/-- … -/` → `theorem …` の順にする。

## 逆典の逐語引用は**記憶で書かない**

**失敗形**: docstring の `原文 (GenEll p.15):` の行を記憶で書くと
`check.mjs` が

    引用照合: 逐語が GenEll 物理 p.15 に見つからない(layout で 31/35 文字まで一致)
    次に来るはず: ")Let"

で NG を出す。★`lake build` は通るので**ビルドだけでは気づかない**。

**直し方**: 同じ項目を引いている**既存のファイルから丸ごと写す**。

```bash
grep -rn "^> Lemma 3.2" lean/ABC3/ --include=*.lean
```

★コミット前に必ず `node tools/check.mjs` を通すこと。

## `addOrderOf_eq_one_iff` は `AddMonoid.` 付き（第 1201）

**失敗形**: `exact absurd (addOrderOf_eq_one_iff.mp h1) h` →
`Unknown identifier 'addOrderOf_eq_one_iff.mp'`。

**理由**: mathlib の `orderOf_eq_one_iff` の `to_additive` 名は
`AddMonoid.addOrderOf_eq_one_iff` である（名前空間が付く）。

**直し方**: `AddMonoid.addOrderOf_eq_one_iff.mp h1`。

## 部分群の membership 証明に `rw [pow_one] at h` は通らない（第 1201）

**失敗形**: `h := (f 1).2`（型は `↑(f 1) ∈ torsionPoints W (l ^ 1)`）に
`rw [pow_one] at h` → `motive is not type correct`
（`l ^ 1` が `f 1` の**型の中**にも現れるため）。

**直し方**: 求める形を `have h : (l ^ 1) • P = 0 := (f 1).2` と
**先に書いて**から `simpa using h`。指数の正規化は `simp` に任せる。

## 在庫を引かずに新規ファイルを作ると名前衝突でビルドが落ちる（第 1201）

**失敗形**: `Found/GaloisRep/TateProjOne.lean` に `tateProj` と
`tateProj_galTate` を新規に書いたら、`ABC3.Found` の import で
`environment already contains 'ABC3.Found.GaloisRep.tateProj_galTate'
from ABC3.Found.GaloisRep.TateWiring`。

**理由**: `tateProj` は `Interface/GaloisRep/Torsion.lean` に
（層 `n` つきの一般形で）、`tateProj_galTate` は
`Found/GaloisRep/TateWiring.lean` に**すでにあった**。

**直し方**: 新しい補題を書く前に必ず
`node tools/decl-index.mjs` → `grep .cache/decl-index.txt <名前>`。
CLAUDE.md 在庫の規則そのものである。★同じ穴は第 1191 でも落ちた。

## 同名の `.src` が 2 つあると check.mjs は無関係なファイルで NG を出す（第 1201）

**失敗形**: `tateProj.src` を重複して定義したら、
`NG lean\ABC3\Found\SemiAnbd\TemperedGroup.lean:129
 G1 tateProj.src の中身を読めなかった` という**まったく別のファイル**の NG。

**直し方**: 重複を消せば直る。★NG の行番号を信じて別ファイルを疑わないこと。

## push の前に `check.mjs` の末尾を読む（第 1184・第 1201）

**失敗形**: `NG 2 件` のまま commit / push した。二度目である。
**直し方**: `node tools/check.mjs 2>&1 | tail -3` が `PASS` を出すまで commit しない。

## `tateModule W l` と `limTors W.toAffine.Point l` は定義から同じ（第 1203）

**測ったこと**: `Interface/GaloisRep/Torsion.lean` の `tateModule W l`
（`torsionPoints W (l^m)` の逆極限）と
`Found/GaloisRep/TateLimit.lean` の `limTors A l`
（`(nsmulHom A (l^m)).ker` の逆極限）は、
membership が `mem_ker_nsmulHom : x ∈ (nsmulHom A m).ker ↔ m • x = 0` で
`Iff.rfl` なので**定義から等しい**。

**使い方**: `limTors` で書いた在庫の補題が `tateModule` にそのまま当たる。
実例: `exact exists_smul_of_proj_zero l n f h`（第 1203、変換なしで通った）。
★抽象側で証明した補題を探すときは `limTors` でも grep すること。

## `ℤ_[l]` に値を取る `def` は `noncomputable`（第 1204）

**失敗形**: `def redVec (w : Fin 2 → ℤ_[l]) : Fin 2 → ZMod l := fun i => PadicInt.toZMod (w i)`
→ `failed to compile definition, consider marking it as 'noncomputable'
because it depends on 'PadicInt.instCommRing'`。

**直し方**: `noncomputable def`。★`ℤ_[l]` が絡む `def` は既定で付けておく。

## `rw [e.apply_symm_apply]` は最初の項でパターンを固定する（第 1205）

**失敗形**: 目標が
`... e (e.symm w) ... = l • e (e.symm u)` のとき
`rw [e.apply_symm_apply]` は `?x := w` で固定され、
右辺の `e (e.symm u)` は**残る**。

**直し方**: 当てたい項を明示する —— `rw [e.apply_symm_apply u]`。
★`rw` は「最初の一致でメタ変数を決め、その具体形をすべて書き換える」規則である。

## `↥M`（`IntermediateField`）の `DecidableEq` は `Subtype` が横取りする（第 1207）

**失敗形**: `M : IntermediateField L L̄` として
`rw [rhPoint_nsmul (algebraMap ↥M L̄) (W.baseChange ↥M) Q' l]` が
「パターンが見つからない」で落ちる。`pp.explicit` で見ると、
在庫側は `Affine.Point.instAddCommGroup … (fun a b ↦ Classical.propDecidable (a = b))`、
こちら側は `Subtype.instDecidableEq …` 由来で、**群構造が別物**になっている。

**理由**: `open scoped Classical` の `Classical.propDecidable` は優先度 low。
`↥M` は `Subtype` なので `Subtype.instDecidableEq`（通常優先度）が先に当たる。

**直し方**: 命題の**中**で `letI` を入れる（`M` が ∃ で束縛されていても書ける）:

```
∃ M : IntermediateField L Lbar, FiniteDimensional L M ∧
  letI : DecidableEq (M : Type) := fun a b => Classical.propDecidable (a = b)
  ∃ Q' : (W.baseChange M).toAffine.Point, …
```

証明の側でも `obtain` 直後に同じ `letI` を置く。★第 1151 の
`HasLCyclicVelu` と同じ穴である（`E.fld` が `↥E.K` だった）。

**診断のしかた**: `set_option pp.explicit true in` を付けた `example` で
同じ `rw` を書き、パターンと目標の**インスタンス項**を並べて見る。

★**別の顔（2026-09-06、第 1444）**: `rw` ではなく `exact` で在庫を当てると、
「パターンが無い」ではなく
`(deterministic) timeout at isDefEq, maximum number of heartbeats (200000)` になる。
`Subtype.instDecidableEq` 由来の `Finset.image` と `Classical.propDecidable` 由来の
`Finset.image` を kernel が突き合わせようとして止まるためで、
**`maxHeartbeats` を上げても本質は直らない**（上の `letI` を入れると 0.5 秒で通る）。
実例: `veluQuotientFull_baseChange (algebraMap K ↥M) E E' hQ hE'`。

## `Ideal.ramificationIdx_le_finrank` は `S K L` が明示引数（第 1209）

**失敗形**: `Ideal.ramificationIdx_le_finrank P.asIdeal` →
`failed to synthesize CommRing ↥P.asIdeal`（`P.asIdeal` が `S` の位置に入った）。

**直し方**: `Ideal.ramificationIdx_le_finrank (𝓞 L') L L' P.asIdeal`。
★`variable (S)` と `variable (K L : Type*)` が効いている。

## `NoZeroSMulDivisors` は `Module.IsTorsionFree` とは別クラス（第 1209）

**測ったこと**: mathlib は `NoZeroSMulDivisors.iff_algebraMap_injective` を
`Module.isTorsionFree_iff_algebraMap_injective` の deprecated alias にしたが、
**`NoZeroSMulDivisors` 自体は残っており別クラス**である
（`.2` の結果は `Module.IsTorsionFree` になって型が合わない）。

**直し方**: クラスの定義（`eq_zero_or_eq_zero_of_smul_eq_zero`）から直接作る:
`⟨fun {c x} hcx => …⟩`（`Algebra.smul_def` ＋ `mul_eq_zero` ＋ 単射性）。

## `io.open(p,'w')` は書き込み前にファイルを空にする（第 1216）

**失敗形**: Python で `io.open(p,'w',encoding='utf-8').write(s)` を使い、
`s` にサロゲート（`'\ud835\udd3d'`）が混じっていて `UnicodeEncodeError` になった。
**ファイルは 0 バイトに truncate されていた**（`open` が先に切り詰めるため）。

**直し方**: 必ず一時ファイル経由にする——
`io.open(p+'.tmp','w',encoding='utf-8').write(s)` → `os.replace(p+'.tmp', p)`。
★これなら書き込みが失敗しても元のファイルは無傷である。
☆復旧は `git checkout -- <file>`（コミット済みなら）。

★サロゲート自体の対処は「Python surrogate escapes」の項——`'\U0001D53D'` と書く。

## `map_adjoin` は `AlgHom`、`map_top` は `Algebra` の名前空間（第 1217）

**失敗形**: `Algebra.map_adjoin` / `Subalgebra.map_top` → `Unknown constant`。

**直し方**: `AlgHom.map_adjoin`（`Mathlib/Algebra/Algebra/Subalgebra/Lattice.lean:865`、
`namespace AlgHom` の中）と `Algebra.map_top`（同 263、`namespace Algebra` の中）。
★`grep -n "^namespace\|^end "` で**行番号より前の namespace** を見て決めること。

## `IsScalarTower (𝓞 L) L M` は自動では出ないが `rfl` で出る（第 1222）

**失敗形**: `M : IntermediateField L L̄` が `L` 上有限次でも
`haveI : IsScalarTower (𝓞 L) L M := inferInstance` が
`failed to synthesize`。`IsScalarTower (𝓞 L) (𝓞 M) M` も同じ。

**測ったこと（2026-09-02）**: 自動で出るのは
`NumberField M`（`NumberField.of_module_finite L M`）・`IsScalarTower ℚ L M`・
`Algebra (𝓞 L) M`・`Module.Finite (𝓞 L) (𝓞 M)`・`Algebra.IsIntegral (𝓞 L) (𝓞 M)`。

**直し方**: 2 つとも
`IsScalarTower.of_algebraMap_eq fun _ => rfl` で出る
——`Algebra (𝓞 L) M` は `𝓞 L ⊆ L → M` の制限だからである。
★**数学の穴ではなくインスタンス探索の経路の問題**。
☆`Found/GaloisRep/TowerInstances.lean` に補題として置いた。

## `RingHom.map_det` の右辺は `mapMatrix` の形（第 1230）

**失敗形**: `RingHom.map_det f M : f M.det = (f.mapMatrix M).det` なので、
`(M.map f).det` を期待して `rw [hmap]` すると
「パターンが見つからない」。

**直し方**: 間に `RingHom.mapMatrix_apply` を挟む
——`rw [hz, RingHom.mapMatrix_apply, hmap] at hdet`。

## `simp` に `redVec_nsmul` を入れると `l •` が `↑l *` に化ける（第 1234）

**失敗形**: `simp only [..., redVec_nsmul, ...] at hkey` の後、
`l • redVec (e u)` が `↑l * redVec (e u)` になり、
`ZMod.natCast_self` を simp に足しても当たらない。

**直し方**: `redVec_nsmul` を simp 集合から外し、
`redVec_nsmul_self : redVec l (l • w) = 0` を直接当てる。
★仮説側の `2 • y` は `y + y` の形で書いておくと `nsmul` が出てこない。

## 仮説の中で「その項の `IsElliptic`」を使うには `∀ [_inst : …]`（第 1248）

**失敗形**: 仮説を
`∀ S, S.card + 1 = l → (velu S).IsElliptic → (… jExp p (velu S) …)`
と書くと、`jExp` が要求する `(velu S).IsElliptic` を
**合成できない**（`→` で受けた命題はインスタンスにならない）。

**直し方**: インスタンス束縛にする——
`∀ S, S.card + 1 = l → ∀ [_inst : (velu S).IsElliptic], (… jExp p (velu S) …)`。
★呼ぶ側は `hrel S hcard hss`（インスタンスは自動で埋まる）。

## `PowerSeries.C` の数値は `map_ofNat` で素の数値に直す（第 1254）

**失敗形**: `ring` が `… * 48 = … * PowerSeries.C 48` を閉じられない
（`C 48` と数値 `48` は `ring` では同一視されない）。

**直し方**: `simp only [map_neg, map_ofNat]` を先に当てて
`PowerSeries.C (-5)` → `-(5 : PowerSeries ℤ)` の形にしてから `ring`。
☆在庫の補題が `C 12 * x` の形で述べられているときは
`simpa [map_ofNat] using h` で数値の形に直す。

## `lake build ABC3` は低層のファイルを触ると 10 分を超える（第 1254）

**失敗形**: `TateSeries.lean`（低層）を編集した後の `lake build ABC3` が
Bash ツールの上限 10 分で打ち切られた（`timeout` に 30 分を渡しても
上限は 600000 ms）。

**直し方**: ビルドだけを単独のコマンドで走らせ、
打ち切られたら**そのまま再実行**する（インクリメンタルに続きから進む）。
☆ゲートと commit は別のコマンドに分ける。

## 在庫を引くときは「探す語」を疑う（第 1256）

**失敗形**: `h4`・`h6`（Tate 曲線の `c₄`・`c₆` の同種関係）を
「eisenstein」で grep して「未着手」と判断した。
実際は `Skeleton/GenEll/TateIsogeny.lean` に
**`c4_velu_tate`・`c6_velu_tate` として証明済み**（sorry 0）だった。

**直し方**: 在庫を引くときは**結論の形**（`c₄`, `c₆`, `veluCurve`, `tateCurveAt`）でも
grep する。☆`node tools/decl-index.mjs` の索引を
「定理が何を言うか」の語で引くこと。
★`grep -c sorry <file>` でそのファイルが完成しているかを先に見る。

## `le_or_lt` は無い(2026-09-02、第 1269)

**失敗形**: `rcases le_or_lt 0 (jExp p W) with hnn | hneg` が
`Unknown identifier 'le_or_lt'` で落ちる(続いて `rcases` が「帰納型でない」と言う)。

**直し方**: `by_cases h : a < b` を使う。分岐の順（`positive` が先）が入れ替わるので、
`·` の中身も入れ替えること。`le_or_gt` / `lt_or_ge` は在るが、
`by_cases` なら名前を覚えなくてよい。

★**2026-09-08 追記（第 1073、`WildJumpChain`）**: 同じ罠にまた落ちた。
★上の字面 `Unknown identifier 'le_or_lt'` は**引用符が違うので grep でも機械照合でも当たらない**。
★現行の Lean が印字する**逐語**はこちらである:

```
error(lean.unknownIdentifier): Unknown identifier `le_or_lt`
Tactic `rcases` failed: `x✝ : ?m.54` is not an inductive datatype
```

★2 行目が必ず伴う（`rcases` の引数が elaborate できないため）。
★`rcases le_or_lt a b with h | h` → `rcases le_or_gt a b with h | h` で**そのまま通る**
（分岐の順は `≤` が先のままなので `·` の中身を入れ替えなくてよい。★`by_cases` より安い）。

## `omit ... in` は docstring の**前**に置く(2026-09-02、第 1271)

**失敗形**: `/-- doc -/` の直後に `omit [Inst] in` を挟むと
`unexpected token 'omit'; expected 'lemma'` になる。docstring は宣言に直結するため。

**直し方**: `omit`／`set_option ... in` は docstring より前の行に置く。

## ゲートの読み方——`grep -A n | head -m` は error を隠す(2026-09-02、第 1271)

**失敗形**: `lake build <target> 2>&1 | grep -E "error|warning" -A12 | head -45` で
「error 無し」と判断したが、他ファイルの warning が 45 行を食い尽くしていて
**自分のファイルの error が表示されていなかった**。そのまま commit した。

**直し方**: 判定は `grep -E "^error" -A 20 | head -40`（行頭アンカー）で行い、
かつ**必ず `lake build ABC3` の tail を読む**。warning と error を同じ grep に入れない。

## `Point.map` の値の型は `W⁄K`(2026-09-02、第 1273)

**失敗形**: `Point.map (S := R) σA P` を
`(W.map (algebraMap R K)).toAffine.Point` の元と足そうとすると
`failed to synthesize HAdd ((W)⁄K).Point ((W.map (algebraMap R K)).toAffine.Point) ?m`。
`W⁄K = W.baseChange K = W.map (algebraMap R K)` は defeq だが構文的に違う。

**直し方**: `noncomputable def tatePointMap ... : A →+ A := Point.map (S := R) σA` のように
**目的の型で一度名前を付ける**。以降はその名前を使えば型は揃う。

## `TateSetup R I K` は `K = Frac(R)` を強制する(2026-09-02、第 1275)

**失敗形**: 同変性 `tatePhi_pointMap` に `σA : K →ₐ[R] K` を渡そうとしたが、
**恒等射しか作れない**。`hmem0`(v ≥ 0 の元は R から来る)と
「x か x⁻¹ の一方は v ≥ 0」から `K = Frac(R の像)` になり、
R-代数準同型は逆元を保つので K 全体で恒等になる(`tateSetup_algHom_eq_id`)。

**直し方**: 環を動かす自己同型が要るときは `tatePhi_map`(σR : R →+* R と σK : K →+* K の対)
を使う。`Point.map` を使いたければ、曲線を σ で固定される**部分環 R₀ の上**に置く。

**教訓**: 仮説の集まりが**その型の元を 1 つに決めてしまう**ことがある。
「この形の σ は本当に非自明に取れるか」を、使う前に 1 度証明して確かめる。

## 在庫を引くときは `head -5` で切らない(2026-09-02、第 1284)

**失敗形**: `grep -rn "pointCoords_rhPoint" ... | head -5` で「無い」と判断して書いたら
`has already been declared`。**5 件の別名(`image_pointCoords_rhPoint_nsmul` 等)が
`head` を食い尽くしていた**。

**直し方**: 宣言の有無を見るときは `grep -rn "theorem <name>\b"` のように
**行頭の宣言キーワードごと**当てる。`head` を付けるなら `-20` 以上にする。

## Point の `+` は DecidableEq のインスタンス経路で割れる(2026-09-02、第 1285)

**失敗形**: `[DecidableEq F]` を variable に置いた file で
`rhPoint_add f W P Q` を `exact` しても
`Point.instAdd` と `Point.instAddZeroClass.toAdd` の不一致で落ちる。
mathlib の `Affine.Point.add` は `x₁ = x₂` の場合分けに `DecidableEq F` を使うので、
**どの DecidableEq を使ったかで `+` が別の項になる**。

**直し方**: 在庫の補題が `open scoped Classical in` の下にあるなら、
**こちらも `open scoped Classical` にして `[DecidableEq F]` を宣言しない**。
どちらを使っているかは在庫の補題の直前の行を見ればわかる。

## 単体ビルドが通っても `lake build ABC3` は落ちる(2026-09-02、第 1288)

**失敗形**: `lake build ABC3.Found.GenEll.GalActVc` が通ったので commit したが、
`lake build ABC3` は
`environment already contains 'ABC3.Found.GenEll.vcPoint_ne_zero.src'` で落ちた。
**同名の宣言が、自分のファイルが import していない別ファイルにあった**ためで、
単体ビルドでは衝突が見えない。

**直し方**: commit 前に**必ず `lake build ABC3` の tail を読む**。
新しい名前を付ける前に `grep -rn "theorem <name>\b" --include=*.lean lean/ABC3`
（`.src` も含めて）当てる。衝突したら在庫の方を import して使う。

## `ZMod (l^1)` と `ZMod l` はイデアルで繋ぐ(2026-09-02、第 1295)

**失敗形**: `rw [pow_one] at h` で `h : (toZModPow 1 D).val < l ^ 1` を書き換えようとすると
`motive is not type correct`——`l^1` が `ZMod (l^1)` の型指数に現れるため。

**直し方**: 型指数には触らず、
(1) 数の不等式は `have hpow : l ^ 1 = l := pow_one l` を足して `omega`、
(2) 環準同型の同一視は**イデアルの言葉**で回す:
`ker (toZModPow 1) = span {l^1}`、`span {l^1} = span {l}`(要素の書き換えなので安全)、
`maximalIdeal = ker toZMod`。

## 「一般形が在庫にある」を先に疑う(2026-09-02、第 1325)

**失敗形**: `vAdd_c4_variableChange`(u の付値 0 なら c₄ の付値は不変)を書いたら
`has already been declared`。在庫（`Found/GaloisRep/LocalHeightDelta.lean:64`）には
**より一般の形**（`vAdd (c₄(C•W)) = vAdd (c₄ W) − 4·vAdd C.u`）が既にあった。

**直し方**: 特別な場合を書く前に、**同じ名前と「一般形」の両方**を grep する。
`vAdd`・`valAdd`・`c₄` のような基本量の変換則は、たいてい一般形で入っている。

## 台帳 `.mjs` の中で Lean 名の `'`（プライム）がクォートを閉じる

`isElliptic_latticeCurve'` のようにプライム付きの宣言名をシングルクォート文字列に入れると `SyntaxError: Unexpected identifier` になる。**プライムを含む行だけダブルクォート**にする。
☆`node tools/_ledger-NNNN.mjs` の出力（`ledger NNNN written`）を必ず目で確認すること——コミットは通ってしまうので気づかない。

## 在庫は「宣言名」ではなく「概念語」で grep する

`sum_veluB_nsmul` のような**これから付ける名前**で `.cache/decl-index.txt` を引いても当たらない。☆探すべきは概念語（`negY`・`stable`・`vcPoint`・`variableChange`）である。
★第 1334 では `veluQuotientFull_variableChange` の 2 仮説を潰す補題を書き上げてから**同名宣言の衝突**（`environment already contains ...`）で在庫の存在に気づいた（`VeluPointSet.lean` 第 949・`VeluImage.lean` 第 912 に全部あった）。
☆`lake build` の衝突エラーは在庫検索の最後の安全網であり、最初の網ではない。

## `DecidableEq` のインスタンス違いは `subst` で揃える

具体の体（`↑(K : IntermediateField ℚ ℂ)` など）では `DecidableEq` が `fun a b => a.instDecidableEq b` に解決され、
`open scoped Classical` を置いた**変数の体**の定理は `fun a b => Classical.propDecidable (a = b)` を焼き込む。
★両者は defeq でなく、`Point` の `+` ・ `addOrderOf` ・ `Finset.image` が全部ずれる。

☆**直し方**——一般の定理の側に `[inst : DecidableEq L]` を**明示の束縛として付け**、
証明の先頭で

```lean
have hinst : inst = fun a b => Classical.propDecidable (a = b) := by
  funext a b
  exact Subsingleton.elim _ _
subst hinst
```

とする。★`inst` は局所変数なので `subst` が通り、以降は古典的な在庫がそのまま使える。
☆呼ぶ側（具体の体）では合成された具体インスタンスがそのまま入るので何も起きない。
（第 1338-1339 で実際にこれで抜けた。）

## `letI` を**定理の主張の中**に置くのは脆い

`∀ (hell : …) (hss : …), …` の前に `letI : DecidableEq F := …` を置いても、
束縛ごとにインスタンスが再合成されて**途中で食い違う**ことがある
（`intro` で zeta 展開されるタイミングの問題）。
☆**推奨**は別の手——一般の補題の側に `[inst : DecidableEq L]` を明示の束縛で付け、
証明の先頭で `Subsingleton.elim` ＋ `subst` で古典的なものに揃える
（上の項を見よ）。★呼ぶ側が具体でも古典的でもそのまま通る。
☆混在（F は古典的、K は具体）の場合は、証明側で `letI` を先に置いて
合成を古典的に寄せる（`exists_ext_point_of_stable` 第 1346 はこれで抜けた）。

## 構造体を返す `def` の射影ではインスタンスが見つからない

`SSCurve.ext` のような `def` が返す構造体の射影 `(E.ext M₀).W` に対しては、
`W.IsElliptic` のインスタンス合成が失敗する（`def` が reducible でないため）。
☆**直し方**——中身の形 `E.W.baseChange ↑(extField E M₀)` で書く。
★`haveI` で手で入れても、別の場所で再合成されると同じことが起きる。

## `AdjoinRoot` の塔は**自前で合成しない**

`AdjoinRoot f`（`f : R[X]`）には mathlib が

* `instance [CommSemiring S] [Algebra S R] : Algebra S (AdjoinRoot f)`
* `instance [IsScalarTower R₁ R₂ R] : IsScalarTower R₁ R₂ (AdjoinRoot f)`

を持っているので、`S → R → AdjoinRoot f` の塔は**そのまま降りてくる**。
☆ここで `(algebraMap R (AdjoinRoot f)).comp (algebraMap S R) |>.toAlgebra` を
自前で `instance` にすると、`SMul` が `AdjoinRoot.instSMulAdjoinRoot` と食い違い

```
Type mismatch … (instAlgebraLocInt p hl).toSMul
  but is expected to have type … AdjoinRoot.instSMulAdjoinRoot …
```

になる。★**まず何も定義せずに `lake build` して、足りないものだけ足す**こと
（第 1377 でこれに落ちた——結局 4 つの自前インスタンスは全部不要だった）。

## 定理の**主張**が要求するインスタンスは `haveI` では入らない

`theorem foo : … IsLocalRing.maximalIdeal C …` のように主張の中で
`[IsLocalRing C]` を使うなら、証明の中の `haveI` では遅い。
☆**直し方**——先に `instance : IsDiscreteValuationRing C := …` を宣言してしまう
（`IsLocalRing` はその射影で降りる）。★第 1377 でこれに落ちた。

## `IsScalarTower.of_algebraMap_eq'` は引数から型が決まらないことがある

`of_algebraMap_eq' rfl` は `IsScalarTower ?m ?m ?m` になって詰まる。
☆`IsScalarTower.of_algebraMap_eq (R := …) (S := …) (A := …) (fun _ => rfl)` と
明示するか、そもそも上の項のとおり**自前で作らない**。

## `Basis` は `Module.Basis` に改名されている

`have b : Basis ι R M := …` は `Unknown identifier `Basis`` になる。
☆`Module.Basis` と書くか `open Module` する（第 1366）。

## Lean の識別子に `′`(U+2032) は使えない

`K′` は `expected token` になる。☆ASCII の `'`（`K'`）を使う。
ドキュメンテーションコメントの中では `′` でよい（第 1374）。

## 新しいファイルを作る前に**同名が無いか見る**

`Write` は既存ファイルを黙って上書きする。
★第 1383 で `Found/GenEll/SplitDichotomy.lean`（第 982 の既存ファイル）を
上書きし、`isCharNeTwoNF_integralModel` が消えて
`SplitAtCompletion.lean` が壊れた。

☆**見分け方**——`Write` の返事が
『has been updated successfully』（上書き）か
『created successfully』（新規）かを見る。
★さらに `Found.lean` への import 追加を
`if '<名>' not in s` でガードしていると、
**既存の同名 import があるせいでスキップされ**、
上書きしたファイルがそのままビルドされて気づきにくい。

☆**直し方**——`git checkout -- <path>` で戻し、
新しい内容は別名（例: `AlphaSplitDichotomy.lean`）に置く。
★**予防**——`ls lean/ABC3/**/<名>.lean` か `git ls-files | grep <名>` を先に見る。

## `.field` を独立行に置くと関数適用に化ける

```lean
def foo : Bar :=
  (baz arg)
    .someProjection   -- ✗ 次の行に置くと `baz arg` への「関数適用」として構文解析される
```

`error: Function expected at baz arg but this term has type Bar` になる
(`lean_check` の対話セッションでは 1 行で書いていたため気づかず、
`lake build` で初めて壊れて見つかった——第 1467)。
★★**`.field`/`.method` は直前の式と同じ行に置くか、式全体を1行にする。**
複数行にしたいなら `Subgroup.topologicalClosure (Subgroup.normalClosure ...)` のように
**外側の関数を先頭に出す**形へ書き換える。

## `structure ... where` の中で `Type*` を使うと以降のフィールドが壊れる

```lean
structure Foo where
  Pic : Type*                     -- ✗
  picAddCommGroup : AddCommGroup Pic
```

`error: Unknown identifier picAddCommGroup`(`autoImplicit` が `false` なので
暗黙変数として拾えない、という副次エラーつき)。★★**`Type*` の auto-bound
universe が構造体フィールドの構文解析を巻き込んで壊す**(2026-09-04、
`ABC3.Interface.LocProP.EtaleSetup` で実測。最小 5 行で再現)。
★**対処: `universe u` を明示して `Type u` と書く。** `Type*` は構造体の外
(トップレベルの `def`/`theorem` の引数)では問題なく使えている。

## `structure` の `[instance field]` は外部で自動解決されない

```lean
structure Setup where
  X : Type u
  [xGroup : Group X]
  ...
theorem foo (E : Setup) (a : E.X) : ... := ...  -- ✗ E.X の Group インスタンスが見えない
```

`failed to synthesize instance of type class Group E.X` になる。構造体の中では
前のフィールドが後のフィールドの型検査に使えるが、**構造体の外で `E : Setup` を
明示引数に取る新しい宣言を書くと、`E.xGroup` は自動で instance キャッシュに入らない**。
★★**対処**: 型注釈を書かず `def foo := @Setup.foo` の形にする(型は Lean が
フィールドの元の型から推論するので instance 解決が要らない)。`theorem` は
型注釈必須なのでこの手が使えない——`def` を使うこと(第 1469)。

## 19. `lean_start` は存在しない import パスでも「成功」を返す(2026-09-04)

`mcp__abc3-lean__lean_start` に `Mathlib.Data.Rat.Basic`・
`Mathlib.GroupTheory.Subgroup.Pointwise` のような**存在しないモジュール名**を
渡すと、「起動して import を読み込んだ」と成功扱いのメッセージが返るが、
**実際には import が 1 つも読み込まれず、後続の `lean_check` が
`unknown namespace`/`unknown identifier` で全滅する**(自分自身の宣言・
`ABC3.Meta` すら見えなくなる)。

★★**見分け方**: 正常時は import 数に応じて数秒(温まっていれば 2〜8 秒程度)
かかるが、**壊れているときは 2〜3 秒で返ってくる**——速すぎる成功は疑うこと。
確実なのは、`lean_start` 直後に `#check` で ABC3 側の既知の宣言を 1 つ
引いてみること(0.01 秒で判定できる)。

**How to apply**: import パスを新規に足すときは、まず 1 個ずつ足して速度を見る
(怪しい 1 個を孤立させる)か、`lean_start` の直後に必ず軽い `#check` で
実際に読めているか確かめてから本題のコードを書く。存在するパスは
`node tools/decl-index.mjs --mathlib` が生成する `.cache/mathlib-index.txt` の
各行末尾のファイルパス(例: `GroupTheory/Commensurable.lean`)から
`Mathlib.` + パスの `/` を `.` に、`.lean` を外した形で機械的に作れる。

## 20. 丸括弧タプル `(a, b, ...)` は `Type` 専用——`Prop` を混ぜると `Prod.mk` 型不一致(2026-09-04)

`def foo := (a, b, c)` の丸括弧記法は常に `Prod.mk`(`Type u`/`Type v` 専用)
に脱糖される。`a`・`b`・`c` の中に `Prop`(`Function.Injective f` の証明項等)
が1つでも混ざると、「sort `Prop` だが `Type` が期待される」型不一致で落ちる
——全要素が `Prop` の場合でも同じく落ちる(`(a,b) : P ∧ Q` は `⟨a,b⟩` の
angle bracket + 期待型の指定が無いと通らない)。

**さらに罠**: 型不一致を直そうとして `def foo (E : S) : A ×' B ×' ... := ⟨...⟩`
のように **`×'`(`PProd`)で明示的に型注釈を書く**と、今度は
`E.field1 : SomeType E.X`(`X` は `S` の中の抽象 `Type` フィールドで、
`[instance]` フィールド `E.XInstGrp` で `AddCommGroup` 等が与えられている)
の `X` に対するインスタンス探索が **注釈された型を独立に再エラボレートする
過程で失敗する**(「配管」#の元祖の亜種——`@Struct.field` を型注釈無しで
使う回避策と同じ根)。

**How to apply**: `Prop` を含む/含みうるタプルは丸括弧もダメ、型注釈付き
`×'`/`∧` もダメ。**`PProd.mk`/`And.intro` を型注釈無しでネストして直接呼ぶ**
(`PProd.mk a <| PProd.mk b <| ... c`)。各 `field` の型は `E` からの
projection として自動的に(インスタンス込みで)決まるので、独立な型注釈も
インスタンス探索も不要になる。実例: `lean/ABC3/Skeleton/Falt1/Section4.lean`
の `theorem_4_1`/`theorem_4_3`/`theorem_4_5`。

## 21. `Over` 圏の base change を手作りしない——`Over.pullback ⋙ Over.map` +
`overPullbackMap` が MorphismProperty の遺伝を一発で与える(2026-09-04)

`f : A ⟶ B`(`Over S` の射)に対し、ある `MorphismProperty` `P`(有限性・
étale 性等)を保つ base change された射(`(-)_K` のような係数拡大の
射側の作用)を構成したいとき、`Limits.pullback.map`(`pullback.lift` 経由)
で手作りすると、`IsFinite`/`Etale` の遺伝を示すのに
`IsPullback.of_right`/`paste_horiz`/`isoIsPullback` 等の pullback 貼り合わせを
自分で組み立てる必要があり、非常に長くなる(1回失敗して撤回した)。

**代わりに mathlib 自身の `CategoryTheory.Over.pullback (f : X ⟶ Y) :
Over Y ⥤ Over X`(base change の関手)と `CategoryTheory.Over.map (f : X ⟶ Y)
: Over X ⥤ Over Y`(構造射との後合成、`.left` を変えない)を合成して使うと、
`CategoryTheory.MorphismProperty.overPullbackMap`
`(f : S' ⟶ S) [P.IsStableUnderBaseChange] {X Y : Over S} (g : X ⟶ Y)
(H : P g.left) : P ((Over.pullback f).map g).left` が
base change 安定性からの遺伝を1行で与える。** `Over.map` 側は
`(Over.map f).map g).left = g.left`(`simp [Over.map]` で示せる)ので、
合成関手 `Over.pullback f ⋙ Over.map f`(対象では `X ↦ X ×_Y S' →(post f)→
Y` のような「base change してまた元の圏に戻す」係数拡大そのもの)の
射側の性質遺伝もこれだけで閉じる。

**How to apply**: `Over` 圏で base change を伴う関手を構成するときは、まず
`Over.pullback`/`Over.map`(と `overPullbackMap`)で書けないか確認してから
`Limits.pullback.map`/`IsPullback` の手作りに進む。実例:
`lean/ABC3/Found/CorrHyp/SchemeFEt.lean` の `ExtF`/`extFEt`
(`AlgebraicGeometry.Etale`/`IsFinite` の base change 安定性を利用)。

## 22. `Cone`/`Cocone` の `naturality` フィールドを直接埋めると `Functor.const`
の配管で詰まる——`NatTrans` を独立に作ってから束ねる(2026-09-04)

`Cone D`(`D : J ⥤ C`)を構造体リテラル `{ pt := ..., π := { app := ...,
naturality := ... } }` で直接組み立てるとき、`naturality` フィールドの型が
`((Functor.const J).obj pt).map f ≫ π.app j' = π.app j ≫ D.map f`
(`Functor.const` で包まれた形、実質 `𝟙 pt ≫ π.app j' = π.app j ≫ D.map f`)
であるため、`pt` が `Limits.pullback`/`Over.mk` 等の「重い」項のとき、
`𝟙 pt`・`π.app j` 等が「`instances` 透明度で型が合わない」で `rw`/`simp`/
`show`/`congr 1` のどれもが詰まることがある——1箇所直しても**別の項で
同じ症状が再発する**(`Over.mk` 単体の既知の罠より根が深い、複合的な配管)。

**How to apply**: `Cone` の `π` を**直接埋めずに**、まず独立した
`NatTrans ((Functor.const J).obj pt) D`(あるいは `D₁ ⟶ D₂` という
関手間の自然変換)として構成する——`Functor.const` の包みが**そこにしか
出てこない**ぶん配管が軽くなり、`show`(展開後の生の等式)+ 個々の補題
(例: `pullback.lift_snd`)で素直に閉じることが多い。それでも `Cone`
リテラルの `naturality` フィールド自体に埋め込む段でまた詰まるなら、
`CategoryTheory.Limits.Cone.postcompose`(`NatTrans (D₁ ⟶ D₂) → Cone D₁ ⥤
Cone D₂`)経由で `Cone` を**自動的に**変換する道を検討する——独立に作った
`NatTrans` を `(Cone.postcompose η).obj s` に渡せば、`Cone` 側の
naturality の再証明が不要になる(`rfl` で確認できることが多い)。
実例: `lean/ABC3/Found/CorrHyp/ExtLimit.lean` の `extDiagram_map_snd`
(独立補題)→ `extDiagramToSpecK`(独立 `NatTrans`、これは閉じた)。

★★★**2026-09-04、頂点側(`π` の定義域が `(Functor.const J).obj pt` 自体)も
含めて完全解決**(`isLimit_extCone`)。残っていた2つの追加の技:

1. **`Cone` を組み立てる `π` フィールド自体は、独立に作った `NatTrans`
   (`extConePi` のように)をそのまま代入する**(`{ pt := ..., π :=
   myNatTrans }`)——構造体リテラルの中で `naturality` を**再度**証明し
   直さない。
2. **`IsLimit.mk` の `uniq` で渡される仮定
   `hm : ∀ j, m ≫ t.π.app j = s.π.app j` は、`rw`/`simp` で直接使おうとすると
   同じ配管に当たる**——`have hm' : ∀ R, m ≫ (myNatTrans).app R = s.π.app R
   := hm` のように**展開後の型を明示して再束縛する**と、`rw` の構文一致
   ではなく `have` 自体の defeq チェック(`set_option backward.isDefEq.
   respectTransparency false` と組み合わせる)を経由するため通る。これは
   「配管の万能薬」と呼べる型の技——`rw`/`simp`/`congr 1`/`show` のどれもが
   `Functor.const` 絡みで詰まったときは、まずこの「`have` で型を明示して
   再束縛」を試す。

実例: `lean/ABC3/Found/CorrHyp/ExtLimit.lean` の `isLimit_extCone`
(`extCone X` が `Ext X` の極限表示であることの完全な証明、§4 `Lemma 4.1`
の構成的降下の核心部品)。

## 23. `abbrev` で定義された `Algebra`/`Instance` は typeclass 探索が
自力で見つけない——`letI := TheAbbrev args` で明示的に呼ぶ(2026-09-04)

`differentIdeal_ne_bot`(mathlib)の暗黙引数 `[Algebra.IsSeparable
(FractionRing A)(FractionRing B)]` を満たそうとして、`Algebra
(FractionRing V)(FractionRing W)` を自分で `RingHom.toAlgebra` で
手作りしたところ、`exact differentIdeal_ne_bot` が
`Algebra.IsSeparable ...` を「見つからない」で落ちた——`haveI hsep :
Algebra.IsSeparable (FractionRing V)(FractionRing W) := by ...` を
直前に置いているのに、である。`@differentIdeal_ne_bot` で引数を
全部埋めて調べると、mathlib 側が要求する `Algebra (FractionRing V)
(FractionRing W)` は自作の instance ではなく
**`FractionRing.liftAlgebra V (FractionRing W)`**(mathlib の
`abbrev`、`Algebra V K → Algebra (FractionRing R) K` を localization
の普遍性から作る既製品)だった——**同じ型でも別の項なので `hsep` が
不一致になり instance search が失敗する**(diamond)。`abbrev` は
`instance` と違い typeclass 探索が自動では見つけない(reducible では
あるが登録されない)ため、**`letI := FractionRing.liftAlgebra ...` で
明示的に呼んで先に `letI`/`haveI` チェーンに乗せる**必要がある。

**How to apply**: mathlib の補題が要求する instance が
`Unknown constant`ではなく「見つからない/型が合わない」で落ちるときは、
まず `@lemma_name` に全引数を渡して**mathlib がどの instance を
選んでいるか**を確認する(エラーメッセージの型注釈に出る)。それが
`abbrev`(`#print` で `abbrev` と出る)なら、自分で `RingHom.toAlgebra`
等を手作りするのではなく**その `abbrev` を `letI` で直接呼ぶ**。
実例: `lean/ABC3/Found/Falt1/Lemma11.lean` の
`falt1_differentIdeal_ne_bot`。

## 24. `class` の `field : ∀ (explicit) [instance] {implicit} ...` で
「明示引数」を渡し忘れると全引数が1つずつズレる(2026-09-04)

`AlgebraicGeometry.Etale.etale_appLE` のシグネチャは `∀ {X Y : Scheme}
(f : X ⟶ Y) [self : Etale f] {U : Y.Opens}, IsAffineOpen U → ∀ {V :
X.Opens}, IsAffineOpen V → ∀ (e : V ≤ f ⁻¹ᵁ U), (f.appLE U V e).hom.Etale`
——`f` が `[self : Etale f]`(instance)より**前の明示引数**になっている。
`Etale.etale_appLE hU hV le_rfl`(`f` を省略、`[Etale f]` を Lean が
勝手に unify してくれると期待)と書くと、`hU` が `f` の位置の
メタ変数を単一化しようとして失敗し、以降 `hV`→`U` の `IsAffineOpen` 引数、
`le_rfl`→`V` の `IsAffineOpen` 引数、という具合に**全部の引数が1つずつ
ズレた場所に入る**——結果、`Type mismatch: Etale.etale_appLE ?m hV ?m
has type ... but is expected to have type Algebra.Etale ...` という、
一見無関係に見える型エラーになる。

**How to apply**: `∀ (explicit) [instance] ...` の並びを持つ `class` の
フィールドを `open` 済みの短い名前(`Etale.etale_appLE` 等)で呼ぶときは、
エラーメッセージが「型が合わない」を返してきたら**まず `#check
@full.name` でシグネチャの引数の並び(どれが明示・インスタンス・暗黙か)
を確認する**。明示引数が `[instance]` より前にあるなら、それを省略せず
必ず明示的に渡す(`Etale.etale_appLE α hU hV le_rfl`)。
実例: `lean/ABC3/Found/CorrHyp/ExtLimit.lean` の `Etale.algebraEtale_appLE`。

## 25. 第22項の追加教訓——型を合わせるための再束縛は `have` ではなく
`set`(型注釈つき)でやる(2026-09-04)

第22項の「`have h' : <明示的に展開した型> := h` で defeq チェックに
迂回する」という技を、`pullback.hom_ext` に渡す2つの射(`m`・`n`、
`Functor.const`/`extDiagram.obj` 越しの非簡約な型を持つ)に適用しようと
したところ、`have hmty : <明示型> := m` という束縛**自体は通る**のに、
直後に `exact hm'`(`hm' : m ≫ ... = ...`、`m` について述べた既存の証明)
を `hmty ≫ ... = ...` の証明として使おうとすると型不一致になった——
`have` は(`Prop` では証明無関係性で問題にならないが)**`Type`(ここでは
`Hom`)の値については中身を消してしまう**(`let`/`set` と違い透明ではない)
ため、`hmty` と `m` が同じ値であることを後続の項が使えない。

**How to apply**: 型注釈つきで再束縛したい対象が `Prop` の証明ではなく
**データ(`Hom` 等の `Type` の項)**のときは、`have x : T := v` ではなく
**`set x : T := v with hx`** を使う(`hx : x = v` が付き、かつ `x` は
`v` に対して透明)。実例: `lean/ABC3/Found/CorrHyp/ExtLimit.lean` の
`extConePi_app_eq`(`m`/`n` を `set ... : <明示型> := ... with h` で
束縛してから `pullback.hom_ext` に渡した)。

## 25. `PUnit`(引数無しの型パラメータ)を補助的に使う定義を `Type*` 化すると、
`PUnit` 自身の宇宙だけが解決不能なメタ変数として残る(2026-09-04)

`MvPolynomial.uniqueAlgEquiv R σ : MvPolynomial σ R ≃ₐ[R] R[X]`(`[Unique σ]`)
を経由して「単変数版」の同型を作る際、`σ := PUnit`(具体の `Type` = 宇宙 0)
で書いている間は問題ないが、囲む定義の `R C` を `Type*`(宇宙変数化)に
した瞬間、`PUnit` の宇宙が**それとは独立の**メタ変数 `?u.48` として残り、
`declaration ... contains universe level metavariables` で失敗する——
`R C` の宇宙をいくら具体化・注釈しても直らない(`PUnit` 自体の宇宙が
別変数だから)。

**How to apply**: `PUnit` を補助的な添字型として使う箇所では、
`(PUnit : Type)`(または必要な宇宙を明示した `PUnit.{0}`)と**書いた
その場で宇宙を固定する**。`σ` を暗黙引数に取る補題(`uniqueAlgEquiv`・
`algebraTensorAlgEquiv` 等)を呼ぶときは `(σ := (PUnit : Type))` で
明示すること。囲む定義の他の型変数(`R C`)を `Type*` にしても、
`PUnit` の宇宙は他と無関係にずっと `Type` に固定されたままでよい
——実際に混在させて(`R C : Type*` かつ `PUnit : Type`)問題なく
`lean_check` が通ることを確認済み。
実例: `lean/ABC3/Found/Falt1/KaehlerAux.lean` の `tensorPolynomialAlgEquiv`。

## 26. `Polynomial.mapRingHom f` を `FunLike` 適用した形と `.map f`(dot記法)
は定義上等しい(`rfl`)のに構文上一致せず `rw`/`simp_rw` が刺さらない(2026-09-04)

`(monomial n a).map f = monomial n (f a)`(`Polynomial.map_monomial`)を
「多項式の多項式」(`Polynomial (Polynomial R)`)の**外側**の階層に適用すると、
`f := Polynomial.mapRingHom φ`(内側の係数環を写す束縛`RingHom`)についての
`f a`(`FunLike`適用、`a : Polynomial R`)が出現する。これは
`a.map φ`(`Polynomial.map`のdot記法)と**定義上完全に等しい**
(`example : ⇑(Polynomial.mapRingHom f) = Polynomial.map f := rfl` が通る)
にもかかわらず、`Polynomial.map_sum`(dot記法`.map`前提)や、dot記法で
書いた別の補題(`key2`等)を`rw`/`simp_rw`で当てようとすると
「instances 透明度で type-correct でない」失敗になる——`FgSubalgebra`
(第22・25項)と同種だが、今回は`FgSubalgebra`の透明度ではなく
**`FunLike`適用 vs dot記法**という別の構文不一致が原因。

**How to apply**: `have hcoe : (⇑(Polynomial.mapRingHom φ) : Polynomial R →
Polynomial S) = Polynomial.map φ := rfl` を明示的に挟んで `rw [hcoe]` して
から先に進む——以後は一貫して`.map`(dot記法)側の補題だけを使う
(`Polynomial.map_sum`・`Polynomial.map_monomial`等、生成的な`map_sum`
ではなく`Polynomial.map_sum`を選ぶ)。実例:
`lean/ABC3/Found/CorrHyp/FieldLimit.lean` の
`exists_fg_subalgebra_tensor_bivariate_finset`。

## 27. `open ... in` / `set_option ... in` はdocstringの**前**に置く——
docstringの後に置くと「unexpected token 'open'; expected 'lemma'」
(2026-09-04)

`/-- docstring -/` の直後に `open X in` や `set_option foo in` を置くと
(`docstring` → `open ... in` → `theorem` という順序)、次の宣言の手前で
パーサが `open`(または `set_option`)を見て「'lemma' を期待している」と
いう紛らわしいエラーを出す——docstringは宣言の**直前**にしか付けられず、
`open ... in`/`set_option ... in` のような修飾コマンドは**その外側**
(docstringより前)に置く必要がある。正しい順序は
`open X in` → `/-- docstring -/` → `theorem ...`。file 内の既存箇所
(`Bivariate_equivMvPolynomial_map` 等)は元々正しい順序だったが、
新規追加時に順序を逆にしてしまい、`lake build` でしか検出されなかった
(`lean_check` は宣言単体を独立に検査するため、この手のファイル内の
前後関係バグは拾わない)。

**How to apply**: `open`/`set_option`付きの宣言を追加するときは必ず
「修飾コマンド→docstring→宣言」の順を確認する。`lean_check`で個々の
宣言が通っても、ファイルへ実際に書き込んだ後は`lake build`(または
該当ファイルだけの`lake build <module>`)で必ず再検査すること——
この種のバグは`lean_check`だけでは検出できない。

## 28. `letI` で導入した `Algebra` インスタンスの下で、既存の等式を
`algebraMap` 形へ変換するときは `▸`/`show` ではなく**ただの `:=`**
(defeq)で通す(2026-09-04)

`letI : Algebra R K := φ.toAlgebra` としてから、別の場所で得た等式
`h : P.map φ = Q`(`φ` を明示形で書いたもの)を`algebraMap R K` を使う
補題(`standardEtalePairPullbackIso` 等)に食わせたいとき、`▸`で
明示形↔`algebraMap`形の間を変換しようとすると(特に`Spec`/`Scheme`の
ような重い型の値に対して)`whnf`がcombinatorial explosionでtimeoutする
——第22・25・26項と同種だが、今回は**変換そのものを避けられる**のが
教訓。`algebraMap R K`は`letI`の下で`φ`と**定義上ぴったり等しい**ので、
`have h' : P.map (algebraMap R K) = Q := h`(型注釈つきの`:=`、`▸`を
一切使わない)だけで通る——`▸`は「型を書き換える」操作としてより重い
defeqチェックを要求するのに対し、`have ... : T := v`は「vがTを満たす
ことを直接検査する」defeqチェックで済み、こちらのほうが軽い場面がある。

**How to apply**: `letI`導入のインスタンスに依存する等式の「形の変換」
が必要になったら、まず`▸`より先に「型注釈つきの`have`/`refine`で
ただ代入できないか」を試す。実例:
`lean/ABC3/Found/CorrHyp/ExtLimit.lean`の`onePieceSchemeIso`
(`hP₀' : P₀.map (algebraMap ...) = Pres.P := hP₀`)。

## 29. opaque な tactic-mode `def`(`AlgEquiv`等)を後から `unfold` して
中の `rw […] at h` 由来の cast を剥がすのは時間対効果が悪い——最初から
「後で使う値」を定義に組み込む(2026-09-04)

`noncomputable def e : A ≃ₐ[R] B := by ... have h := ...; rw [heq] at h;
... exact h.trans ...` のように、`rw [heq] at h`(`heq : x = y`、`x`・`y`
が defeq でない genuine な命題的等式)で中間結果の**型**を書き換えてから
`.trans`/`.comp`する構成をした場合、生成される項には `cast`(`Eq.mpr`/
`Eq.mp` 経由)が埋め込まれる。この `e` を**後から** `unfold e; simp only
[eq_mpr_eq_cast, eq_mp_eq_cast, cast_eq, ...]` で剥がして「`e` を具体的な
元(例: 生成元)に当てた値」を計算しようとすると、外側の(`set`由来の)
cast は`cast_eq`で消せても、`heq`由来の**内側の**cast(関数型`A ≃ₐ[R] B`
自体への cast)は`cast_eq`・`AlgEquiv.trans_apply`・`Subtype.ext_iff`の
どれでも綺麗に剥がせず、`subst`/`generalize`も(`x`が`set`由来のletだと)
「motive is not type correct」で効かない——各試行が**50〜75秒**かかる
(`whnf`が巨大な項を舐める)のに、複数回試みても収束しなかった
(実例: `falt1AdjoinRootEquivIntegralClosure`から`e(root f)=root g`を
取り出そうとした試み、falt1-goal.mdの2026-09-04分に詳細記録)。

**How to apply**: `e (具体的な元)` の値を**後で**知りたいなら、`e`の
構成そのものを変えて「その値が欲しい形にあらかじめ定義する」方を選ぶ
——例えば `w := e (root f)` と**先に定義してから** `e`のstatementの
方を`w`を使って書き直せば、対応する等式は`rfl`で済む(cast を通す
必要が最初から無い)。「まず抽象的な同型を作ってから、後で `unfold`
して具体的な元での挙動を調べる」という順序そのものが罠——`rw […] at h`
で型を書き換える構成をする**前に**、後で必要になる具体的な等式が
何かを見極め、それが`rfl`になるように定義の**順序**を選ぶこと。

## 29. `noncomputable def`の中で`Exists.choose_spec`を分解するときは
`obtain`ではなく`let`+`.1`/`.2`射影を使う——`obtain`の`And.rec`は後で
`unfold`しても簡約されずスタックする(2026-09-04)

`def foo := by obtain ⟨h1, h2⟩ := someProp.choose_spec; ...`のように
`def`(または`noncomputable def`)の中で`obtain`を使って`Exists.choose_
spec`(`Prop`の`And`)を分解すると、生成される項の内部に`And.rec (fun h1
h2 => ...) someProp.choose_spec`という形の**簡約されないパターンマッチ**
が残る——`someProp.choose_spec`は`Classical.choice`経由の**不透明な証明
項**なので、`And.rec`は(具体的な`⟨_,_⟩`構成子に対してしか)β簡約できず、
後で別の場所からこの`def`を`unfold`して`simp`で内部を書き換えようとしても、
`And.rec (fun ... => 巨大な式) ⋯`という形のまま止まってしまい、`simp`が
中の`Iso.trans`等の構造に一切触れられない。

**How to apply**: `def`の中で`Exists.choose_spec`(や`Iff`・`And`一般)を
分解する必要があるときは、`obtain ⟨h1,h2⟩ := ...`ではなく
`let h1 := (...).choose_spec.1`・`let h2 := (...).choose_spec.2`
(フィールド射影)を使う。後で`unfold`+`simp`する予定があるなら**必ず**
こちらを使うこと——射影は`And.rec`のような`rec`を経由しないため、
`unfold`後に`simp`がそのまま構造の中へ入っていける。実例:
`lean/ABC3/Found/CorrHyp/ExtLimit.lean`の`gdT`(`GlueData`の遷移射)。
同じファイルで`set`が`.choose_spec`由来の仮定への`rw`と噛み合わない
場面にも当たった(第25項の亜種)——`.choose_spec`を扱うときは`set`より
先に「`have`を使わず生の式のまま`rw`する」を試すとよい。

## 30. 不透明な`def`で包んだ射を`pullback.fst/snd`等の引数に渡すと、
`isPullback_opens_inf`系のsimp補題が「未使用」のまま一切効かない——
`@[reducible]`を付けて`instances`透明度で展開できるようにする
(2026-09-04)

`gdF (i j : J) : gdV i j ⟶ Z i := (…).ι`のような**不透明な`def`**(`gdV`・
`gdF`とも普通の`noncomputable def`)を`pullback.fst (gdF i j) (gdF i k)`
のように使うと、`IsPullback.isoPullback_hom_fst`・`_hom_snd`・`_inv_fst`・
`_inv_snd`(と`_assoc`版)を`simp only […]`にいくら渡しても**すべて
「未使用」**になる。`rw`で直接試すと理由が分かる:
```
Application type mismatch: The argument
  pullback.snd (gdF …) (gdF …)
has type
  … (pullback (gdF …) (gdF …)) …
but is expected to have type
  … (pullback U.ι V.ι) …
Note: The target expression is not type-correct under the `instances`
transparency level, which may have triggered the failure.
```
`gdF i j`と`U.ι`は`default`透明度では definitionally equal(`gdF`を
展開すれば同じ)だが、simpの書き換え一致判定は`instances`透明度を使う
ため、`gdF`のような**通常の`def`はそこで展開されない**——結果、
`pullback (gdF i j)(gdF i k)`と`pullback U.ι V.ι`が「別の型」として
扱われ、書き換え全体が失敗する。`unfold gdF`をゴールに対して先に
行っても直らない(`unfold`は表面のシンタックスを書き換えるだけで、
`HasPullback`インスタンス経由で決まる`pullback`対象そのものの
「同じ透明度で同じに見えるか」という判定には影響しない)。

**How to apply**: このように**後で`isPullback_opens_inf`系(または他の
`instances`透明度前提のsimp補題)と組み合わせて使うつもりの`def`**
(特に開埋め込み`.ι`をラップするようなもの)は、最初から`@[reducible]`
を付けておく。トイ例で先に確認してから本体に適用するとよい:
```lean
@[reducible] noncomputable def testF {X : Scheme} (U : X.Opens) :
    (U:Scheme) ⟶ X := U.ι
```
これで`isoPullback_hom_snd`等が問題なく`simp`で適用できるようになる。
副作用: `@[reducible]`化すると、他の場所の`unfold gdV; rfl`のような
明示的な展開ステップが不要になる(自動的に`instances`透明度で展開
されるため`rfl`が`rw`の中で自動的に成立する)——ビルドエラー
("no goals")として現れるので、その`unfold …; rfl`を消せばよい。
実例: `lean/ABC3/Found/CorrHyp/ExtLimit.lean`の`gdV`・`gdF`
(コミット`e41d967e`)。

## 戻り値の型が特定の`Algebra`instanceを直接使う`theorem`は型検査の
## 時点で詰まる——`Σ'`で instance ごと束ねて`def`にする(2026-09-04)

**症状**: `theorem foo ... : @Algebra.IsSeparable R A _ _
(FractionRing.liftAlgebra ...) := by ...`のように、**戻り値の型**
(結論)の中で`FractionRing.liftAlgebra`(または他の`[FaithfulSMul ...]`
等の非自明な instance を要求する構成子)を直接使うと、証明本体
(`by ...`)に入る**前**——シグネチャそのものの型検査の時点——で
`FaithfulSMul`・`Fact (Irreducible ...)`等の instance 検索が走る。
シグネチャの仮定だけではこれらを満たせない場合、
`failed to synthesize instance` や(仮に満たせても)`whnf`/`isDefEq`
の timeout(数十秒〜)に化ける。`mcp__abc3-lean__lean_check`の
断片テストでは(仮定を後から`haveI`で足せてしまうので)気付きにくい。

**直し方**: 結論を`Σ' (inst : Algebra R A), @Algebra.IsSeparable R A _ _
inst`のように、**instance そのものを戻り値として束ねる**
(`Prop`ではなく`Type`になるので`theorem`ではなく`noncomputable def`
にする)。こうすると`inst`は単なる束縛変数になり、シグネチャの型検査は
instance 検索を一切要求しない。呼び出し側は`.1`を`letI`で登録すれば
`.2`がその instance に対する証明として使える。
実例: `falt1_hsep_bundled`(`Found/Falt1/KaehlerAux.lean`、
コミット`cd7a95ec`)。

## 31. `unfold`+`simp`で作った`Scheme.basicOpen`絡みの巨大な項に対して
`rw`/`simp`/`conv`はすべて詰まる——`congrArg`+`(Category.assoc _ _
_).symm`のterm-modeで組み立て、部品は独立した`theorem`に先出しする
(2026-09-04)

**症状**: `Scheme`の`isoImage`/`eqToIso`/`basicOpen`を何層も重ねた
ゴール(`unfold <定義>`+`simp only […]`で作る)に対して、さらに別の
事実(`have step := ...`で正しく型検査できる、ゴールと完全に一致する
主張)を`rw [step]`・`simp only […, step]`・`conv => rw […]`のいずれで
差し込もうとしても、判で押したように次のエラーになる:
```
Application type mismatch: The argument
  X.presheaf
has type
  TopCat.Presheaf CommRingCat ↑X.toPresheafedSpace
but is expected to have type
  (TopologicalSpace.Opens ↥X)ᵒᵖ ⥤ CommRingCat
in the application
  X.presheaf.obj
Note: The target expression is not type-correct under the `instances`
transparency level, ...
```
`set_option backward.isDefEq.respectTransparency false`を足しても直ら
ない。`Scheme.basicOpen`自体がmathlib側の定義で内部的に`X.presheaf`
(`TopCat.Presheaf`、`instances`透明度では`Xᵒᵖ⥤C`へ展開されない)に
依存しているため、`rw`/`simp`/`conv`共通の congruence motive構築
(`kabstract`)がこの型を跨げないのが原因と見られる——挿入する事実自体は
`have`単体なら常に正しく型検査できるのに、**ゴールへの適用だけ**が
一貫して失敗する。

**直し方(唯一有効だった方法)**: `rw`を一切使わず、`calc`+`congrArg`+
`Category.assoc`を**termとして**(`by rw […]`ではなく`:=`で直接)
組み立てる。`congrArg f step`は`step : a = b`と関数`f`から`f a = f b`
を直接構成するだけで、`rw`のような「ゴールの中からパターンを探す」
motive探索(`kabstract`)を経由しないため、この壁に当たらない。
再結合(`(f≫g)≫h = f≫(g≫h)`)が必要な箇所も`rw [Category.assoc,…]`
ではなく`(Category.assoc _ _ _).symm`(または`Category.assoc _ _ _`)
を直接`calc`の等式の証明項として使う——`rw [Category.assoc]`単体でも
同じ「`instances`透明度で型が合わない」エラーになることがある。

**もう1つの罠**: 上の`calc`を1つの巨大な`theorem`の中で、必要な事実
(`eqToIso_homOfLE_comm`の適用結果等)を`have step1 := …; have step2 :=
…; …`と内部で構築しながら書くと、それだけで`whnf`のheartbeat上限
(`set_option maxHeartbeats 4000000`でも、400秒の壁時計タイムアウトでも
不足)に達することがある。**各部品を先に独立した`theorem`として証明し
切ってから、本体の`calc`ではその名前を参照するだけにする**と、劇的に
軽くなる(数十秒→1秒未満の部品もある)——閉じた項を`theorem`として
確定させると、以降の`whnf`はそれを不透明な定数として扱えるため、
`have`で毎回インライン展開されるのと違って計算が繰り返されない。

実例: `lean/ABC3/Found/CorrHyp/ExtLimit.lean`の`transitionElemIso_
inv_naturality`(および部品`transitionElemIso_step1`〜`step45`・
`transitionElem_restrict_mul_le`)、コミット`c4172c85`。3ターン
連続でこの壁に当たり続けた末にたどり着いた対処法。

## ヘルパー補題内部の`haveI`は外に漏れない——`AdjoinRoot f`(条件付き
## `Field`)経由の定理を呼ぶ前に自分でも同じ instance を再構築する(2026-09-04)

**症状**: `AdjoinRoot fK`(`fK`の既約性に依存して`Field`になる型)上の
元`w`について、ヘルパー補題`foo`(内部で`haveI : Fact (Irreducible
fK) := …`等を使って`w`の性質を証明する)を呼んで得た事実`hw`を、
**別の**定理`bar w hw ...`(`bar`もまた`AdjoinRoot fK`が体であることを
要求する)に渡すと、`Application type mismatch`(`hw`の型が`bar`の
期待する型と合わない)+ `(deterministic) timeout at whnf` という
紛らわしい形で失敗する。`hw`の型を目視で見比べても完全に一致して
見えるので「instance の diamond だ」と誤診断しやすい
(実際にこのセッションでも一度そう誤診断した)。

**真因**: `foo`の内部で`haveI : Fact (Irreducible fK) := …`等を
した instance は`foo`の**証明の中だけ**で有効で、`foo`の**戻り値の
型**に現れない限り呼び出し側には伝播しない。呼び出し側(`bar`を呼ぶ
その場所)には`Fact (Irreducible fK)`が存在しないので、`AdjoinRoot
fK`がそもそも`Field`だと分からず、`bar`が要求する
`FiniteDimensional`・`Algebra.IsSeparable`等の instance 探索が
(存在しない前提から)非常に高価な/失敗する探索に迷い込む。

**直し方**: `bar`を呼ぶ**その場所でも**、`Fact (Irreducible fK)`・
`FiniteDimensional K (AdjoinRoot fK)`・`Algebra.IsSeparable K
(AdjoinRoot fK)`を(ヘルパー内部と同じ手順で)`haveI`で再構築してから
呼ぶ。これだけで(`isDefEq`timeoutという症状のわりに)一瞬で解決する
——`@`明示引数や`show`での型強制は的外れな対症療法だった。

実例: `falt1_hspan_eq`(`Found/Falt1/KaehlerAux.lean`、コミット
`ea63551e`)——`differentIdeal_eq_span_derivative`/
`conductor_mul_differentIdeal`を呼ぶ前に3つの instance を
再構築して解決した。

## 32. `unfold`を同じ定数の**2つ以上の出現**に同時適用すると、以降の
`rw`/`simp`/`exact`/`refine`/`show`が軒並み極端に重くなる——1つずつ
名前付きの事実にしてから`unfold`を封印する(2026-09-04)

**症状**: `A i j k`・`A j k i`・`A k i j`のように、**同じ`noncomputable
def`(`A`)の異なる引数での出現が複数**ゴールに現れる状態で`unfold A`を
かけると、そのあとどんな手段(`rw`・`simp only`・`exact`・`refine`・
`show`)で式を差し込もうとしても`maxHeartbeats`を200万→400万→2000万
まで上げても完走しない、という壁に当たる。`A`の出現が**1つだけ**の
ゴールに対する`unfold A`は(同じ`A`の定義でも)何の問題も無く軽い
(0.05〜0.5秒)——出現数が2つになるだけで質的に別の壁に変わる。

`#31`(`instances`透明度の壁、`X.presheaf`絡み)とは**別の現象**で、
エラーメッセージも出ない(単に長時間終わらない、または`whnf`/
`isDefEq`のheartbeat上限に達するだけ)。`unfold`自体(ゴールを画面に
表示するだけの段階)はどちらのケースでも一瞬で終わる——重いのは
**その後**の型検査。

**直し方**: `A x y z`という形の式を、`unfold A`を経由せず**名前だけで
参照できる事実**(`A_eq : A x y z = ⟨明示的な右辺⟩`、`unfold A; simp
only […]`だけで証明する——これは`A`の出現が1つだけなので軽い)として
先に確定させる。複数の`A`出現を扱う本体では、`unfold A`を**二度と
使わず**、`rw [A_eq …, A_eq …, A_eq …]`で明示形に置き換えるだけに
する。これだけで劇的に軽くなる(0.05〜0.5秒)。

もう1つ、この直し方と組み合わせて初めて効果が出た教訓: `congrArg`で
式の一部を書き換えるとき、**書き換え対象を2つの合成の"間に挟む"**
(`congrArg (fun x => A ≫ x ≫ B) h`)と、たとえ`A`の出現が1つしかない
文脈でも重くなることがある——`A`・`B`の型が任意の`x`に対して整合する
ことを確認する型検査(ジェネリックな`x`を挟んだ両側の合成可能性の
検証)が高くつくと見られる。**常に「前だけ」(`congrArg (· ≫ K) h`)
か「後ろだけ」(`congrArg (K ≫ ·) h`)の`congrArg`を順番に適用する**
ことで回避できる——`x`を式の"端"にだけ置き、決して真ん中に挟まない。

実例: `lean/ABC3/Found/CorrHyp/ExtLimit.lean`の`gdT'_pair`/
`gdT'_cocycle`(`Scheme.GlueData`の`cocycle`フィールド、コミット
`9fc6c7ad`)。`gdT'_hom_eq`/`gdT'_inv_eq`(`gdT'`の出現1つだけを
`unfold`して確定させた明示形)を用意してから`gdT'_pair`を`rw`だけで
組み立て、`congrArg`は常に「前だけ」/「後ろだけ」で適用した。

## 33. `instances`透明度の壁(`#31`)は型クラス**探索**にも起きる——
`def`(`@[reducible]`でない)越しの射を渡すと`haveI`が効かない、
明示的な型注釈付き`letI`で正規化してから渡す(2026-09-04)

**症状**: `QcqsFEt A B := FEtK A.1 B.1`のような、`@[reducible]`でない
`def`で包んだsubtype(`{f : A.1 ⟶ B.1 // IsFinite f.left ∧ Etale
f.left}`)の要素`c.α`について、`c.α.1.left`(射の中身)を`corrPieceGlueData
X.1 U hU c.C.1.left c.α.1.left`のように`[IsFinite α][Etale α]`を要求
する関数へ**直接**渡すと、たとえ直前に`haveI := c.α.2.1`(または
`haveI : IsFinite c.α.1.left := c.α.2.1`という型注釈付きの形)で局所
instanceを登録しても、
```
failed to synthesize instance of type class
  IsFinite (Over.Hom.left ↑c.α)
```
で失敗する。原因は`#31`と同じ透明度の壁——`c.α.1.left`の「自然な」型
(`c.C.1.left ⟶ (QcqsExt X).1.left`)と、呼び出し先が要求する型
(`c.C.1.left ⟶ (ExtF.obj X.1).left`、`QcqsExt X := ⟨ExtF.obj X.1,…⟩`
の`.1`を`delta`展開しないと一致しない)が`instances`透明度では別物
に見え、型クラス**探索**(`synthInstance`)が局所instanceを見つけられ
ない——`have`の型注釈を書いても、ELABORATION自体は`default`透明度で
成功するので気づきにくい(型注釈の`have`は素通りし、その`have`を
**使う側**の型クラス探索だけが失敗する)。

**直し方**: 呼び出し先が要求する型を**そのまま構文として書いた**
`letI 変数名 : 要求される型 := 元の式`を用意し、以降はその変数名だけ
を使う——`元の式`自身の型注釈ではなく、**呼び出し先の期待する型で
包み直す**のが鍵:
```lean
letI hα : c.C.1.left ⟶ (ExtF.obj X.1).left := c.α.1.left  -- ここで正規化
letI : IsFinite hα := c.α.2.1
letI : Etale hα := c.α.2.2
corrPieceGlueData X.1 U hU c.C.1.left hα  -- hαを使う、c.α.1.leftを直接使わない
```
これで`hα`の**構文上の型**が呼び出し先と完全一致するので、局所instance
が即座に見つかる。

実例: `lean/ABC3/Found/CorrHyp/Instance4.lean`の`corrPieceGlueDataOfCorr`
/`corrPieceGlueDataOfCorr_cover`(`corrHypInstance4`の`Corr`実データを
`ExtLimit.lean`の`corrPieceGlueData`へ接続する配線、コミット`2471cb91`)。

## 34. `Classical.choice`由来の深いnested `Exists.choose`を経由する項は、
**独立に書いた型注釈と照合させる**(`theorem`+`exact`)と`whnf`が
`maxHeartbeats`を100万まで上げても止まらない——型注釈を省略し
`noncomputable def`でLeanに**推論させる**と一瞬で通る(2026-09-05)

**症状**: `Nat.rec`で無限列を組む際、各段で`Classical.choice`
(`Exists.choose`)を5段ほどネストして使う構成(`「1段分の全射性」
定理の返す∃を5個choose_specで剥がしてstructure literalへ詰める`)
から得た値`v := {pt:=hex.choose, ...}`について、「`v`(または`v`を
経由してさらに別のstructureへ詰め替えたもの)の`.pt`/`.hn`/`.hmem`
を使った**独立な型注釈**を`theorem foo : <型> := v.hcompat`のように
書くと、
```
(deterministic) timeout at `whnf`, maximum number of heartbeats
(1000000) has been reached
```
で刺さる。`set_option maxHeartbeats`を100万(既定の5倍)まで上げても
解決しない——単なる「遅い」ではなく、この形の照合自体が実質的に
終わらない。★★紛らわしい点: **個々の射影の等式**(例:
`(psiGenStep K...).pt = (psiGenStepResult K...).pt`)は`rfl`で
0.3秒程度で通る。刺さるのは「複数の射影を組み合わせた独立な型注釈
を書いて、既存の項がそれに一致するとELABORATORに照合させる」局面
だけ——個々の部品はいくら軽くても、**組み合わせて独立に書き直した
型**との照合は別問題として重い、という非対称な挙動。

**直し方**: 型注釈を**省略**する。`theorem foo : <型> := proof`を
やめて`noncomputable def foo := proof`と書き、Leanに`proof`自身の
型を**そのまま推論**させる(`#check`で見ると`v.hcompat`が実際に
証明している`v`自身のフィールドを経由した型になる——見た目は
遠回りだが、独立な型を照合させる工程が丸ごと消える)。ダウン
ストリームで「もっと素直な形の型」が要る場面が来たら、そこで
初めて**個々の射影の`rfl`**(こちらは軽い)で橋渡しする——大きな
組み合わせ型を一度に照合させようとしないのが鍵。

実例: `lean/ABC3/Found/PGC/LubinTateGeneratorSequence.lean`の
`psiGenStep_compat`・`psiGenSeq_compat`(無限compatible列の構成、
コミット`184a7d60`)。

★★続報(2026-09-05、`LubinTateReciprocityLimitCompat.lean`):
**逆に「型注釈を省略すると通らない」場面もある**——`#34`の教訓を
そのまま適用して`reciprocityMapLimitCompat`の型を省略すると、今度は
`Eq.trans key hcongr`という単純な操作にもかかわらず`?m`という未解決
メタ変数が残ったまま`exact`が失敗した(`Classical.choice`由来の罠
とは**別の**、単に「型注釈なしの`def`+大きな`by`ブロック」が型推論に
失敗するという、より平凡な現象)。★対処: **型注釈を省略せず明示的に
書き**、かつ`reciprocityMap`が要求する`FiniteDimensional`インスタンス
2つを`[...]`の明示的な引数として追加する(`.hfd`経由でしか手に入らな
いため)——これで独立に書いた型注釈が今度は問題なく通った。**教訓**:
「型注釈を省略する」も「明示的に書く」もどちらも銀の弾丸ではない——
`Classical.choice`由来の深いnested chooseが**型注釈の側**にある時は
省略、`FiniteDimensional`インスタンスが**項の側**(`.hfd`)からしか
出せない時は明示、と使い分けが要る。迷ったら両方試す。

## 35. 第一同型定理の`≃*`を`mk`の上で計算する自然性——`rw`ではなく`▸`で
構成し直し、橋渡しの補題は**大域`theorem`として**切り出す(2026-09-05)

**目的**: `principalUnitsQuotientEquiv K hπmax n hn (QuotientGroup.mk u)
= unitReductionQuotientMap K π n u`(`(𝒪_K)^×⧸principalUnits(n)≃*
(𝒪_K/π^n)^×`という第一同型定理由来の同型が、`mk`の上では単なる還元
写像そのものとして計算できる、という自然性)を示したい場面。

**罠その1(`rw`と`▸`は同じ結論でも違う項を作る)**: `principalUnits
QuotientEquiv`はもともとtactic-modeの`rw [principalUnits_eq_ker]`で
ゴールの型を書き換えてから`quotientKerEquivOfSurjective`を`exact`して
いた。この構成のまま`unfold`+`generalize_proofs`+`induction`/`cases`
で自然性を攻めると、`Eq.mpr (congrArg (fun S => ...) h) e`という
cast項の**motive**が`rw`の`kabstract`が選んだものになり、独立に書いた
`h ▸ e`の項(TERM-modeの`▸`が推論するmotive)と**形が食い違う**
(前者は単純、後者はTypeclass引数まで巻き込んだ複雑なΠ型になったり
する)。両者は命題として同じ`h`を使っていて**命題としては同一**でも、
`Eq.rec`は`h`が`rfl`に簡約されない限り計算で潰れないため、`exact`の
defeqチェックがこの食い違いを飲み込めない。★直し方: `principalUnits
QuotientEquiv`**自体の定義**を`rw`ではなく`principalUnits_eq_ker K π n
▸ QuotientGroup.quotientKerEquivOfSurjective _ hsurj`という項モードの
`▸`に書き換える(型は完全に同じなので下流に影響ゼロ)。これで定義側と
橋渡し補題側が**同じelaboration経路**を通るようになり、`subst`一発の
議論が届くようになる。

**罠その2(`Subgroup.Normal`をtelescopeの独立引数にすると`▸`の
motiveが余計に太る)**: 橋渡し補題を`{S:Subgroup G}[S.Normal](h:S=φ.ker)`
という形で書くと、`h▸e : G⧸S≃*H`の型注釈のmotive推論が`[S.Normal]`
という**別の局所仮定**まで一緒に汎化してしまい(`S`に依存する仮定は
`▸`が自動的に道連れにする)、`Eq.ndrec (motive:=fun {S}=>[S.Normal]→
S=φ.ker→G⧸S≃*H) ...`という複雑な形になって、実際のゴール(`(𝒪_K)^×`
はアーベル群なので`Normal`は`∀x,x.Normal`という**値に依らない一様な
インスタンス**から来る)と噛み合わない。★直し方: `[Group G]`ではなく
**`[CommGroup G]`で書く**——アーベル群の部分群はすべて自動的に`Normal`
になる(一様なインスタンス)ので、`[S.Normal]`を独立引数として書く
必要が最初から無くなり、`▸`のmotiveが単純なまま保たれる。

**罠その3(ローカルな`have`束縛の依存関数は、位置引数を明示的に全部
与えても前の引数がメタ変数のまま残る)**: 罠1・2を修正した補題を
```lean
have key : ∀ {G H:Type*}[CommGroup G][Group H](φ:G→*H)(hsurj:...)
    {S:Subgroup G}(h:S=φ.ker)(g:G), (h▸...) (QuotientGroup.mk g) = φ g := ...
exact key φ₀ hsurj₀ h₀ u   -- 4引数すべて具体的な項
```
のように**ローカルな`have`**として立てて4引数すべて明示的な項で
適用しても、`h₀`を型検査する段階で`φ`が`?m`のまま残り、「期待型
`Eq.{u+1} ?m (MonoidHom.ker ?φ)`だが実際は`Eq.{1} ...`」という宇宙
不一致で失敗する——`refine key ?_ ?_ h₀ ?_`でも`(φ:=...)`という
named引数でも症状は同じ。単純化した`ℤˣ`上の例でも同様に再現する
(最小再現は本entryのコミットの差分を参照)。★直し方: **`have`を
やめて大域`theorem`として環境に積む**(あるいはファイルに`theorem`
として書く)。同じ4引数の適用が、大域宣言に対してなら一発で通る——
ローカルな依存関数適用のエラボレーション順序に特有の癖で、宣言を
大域化するだけで消える。

実例: `lean/ABC3/Found/PGC/UnitsInverseLimit.lean`の
`quotientKerEquivOfSurjective_cast_apply`(罠1・2・3すべてを踏まえた
大域補題)・`principalUnitsQuotientEquiv_apply_mk`(それを使う自然性
定理)、および`lean/ABC3/Found/PGC/AdjoinIntegers.lean`の
`principalUnitsQuotientEquiv`(`rw`から`▸`への書き換え)。

## 36. 同じ`algebraMap`/`CommRingCat.ofHom`の生の式を複数箇所に書くと
`pullbackSpecIso`等との単一化に失敗する——名前を1回だけ確定させた
独立`def`を経由する(2026-09-05)

**症状**: `Ideal.Quotient`商(`MvPolynomial ... ⧸ I`)を`Algebra`の
コドメインとして使う`algebraMap R (MvPolynomial ... ⧸ I)`を、ゴールの
`letI`チェーンの中で**生の式のまま**2箇所以上(型の`Nonempty (...)`
本体と証明の`refine ⟨(pullbackSpecIso R S T).trans ?_⟩`)に書くと、
```
Application type mismatch: The argument
  algebraMap R S
has type
  @RingHom R S Algebra.TensorProduct.instCommSemiring.toNonAssocSemiring
    (@Semiring.toNonAssocSemiring S (Ideal.Quotient.semiring I))
but is expected to have type
  @RingHom R (?m ...) CommRing.toCommSemiring.toNonAssocSemiring
    (@Semiring.toNonAssocSemiring (?m ...) CommRing.toCommSemiring.toSemiring)
```
のような、**同じ`S`に対する2つの非`defeq`な`Semiring`/`CommRing`
経路**が衝突するエラーになる。`letI hCR : CommRing (...) := inferInstance`
で基底環の`CommRing`を先に固定しても直らない(`S`側の`CommRing`は
別途その場で再導出されるため)——`#31`/`#33`(instances透明度の壁)の
親戚だが、今回は`rw`ではなく`pullbackSpecIso`のような**強く型付けられた
mathlib補題への直接適用**で起きる新しい失敗形。`set S0 := <式> with hS0`
で名前だけ付けても、ゴール側の元の出現が構文的に一致しないため
置換されず(`unfold`直後だと特に起きやすい)、症状は変わらない。

**直し方**: その`algebraMap`を使うスキーム射自体を**独立した`def`として
1回だけ確定させ**(`standardEtalePairSpecMap`が既にこのパターン)、
以降はその`def`の名前だけを参照する——`pullbackSpecIso R S T`を直接
呼ぶ側は、`S`の位置に生の`Ideal.Quotient`式ではなく、その`def`の
`unfold`で出てくる**同じ1回だけ確定させた式**を使う。要するに
「同じ複雑な式を離れた場所で2回書かない・1箇所で`def`にして名前で
参照する」という原則を、`algebraMap`/`CommRingCat.ofHom`の場面にも
適用する。

実例: `lean/ABC3/Found/CorrHyp/ExtLimit.lean`の`descendPieceR_toBase`・
`descendPieceR_reBaseMap`(`algebraMap`をそれぞれ1回だけ確定させた
独立`def`)・`descendPieceR_iso`(それらの名前だけを参照して
`pullbackSpecIso`を適用)、コミット`65fd0a77`。

## 37. 「別の構造体経由の表現」と「直接の表現」が`rfl`で一致する
ときは、`rw`ではなく`congrArg`/`exact`(defeqチェック)で橋渡しする
(2026-09-05)

**症状**: `PsiGenStepResult`(`#34`で導入した、`Classical.choice`由来の
中間構造体)の`.pt`フィールド経由で書かれた既存の定理(例:
`reciprocityMapLimitCompat`)を、`psiGenSeq (m+1)`という**別の(だが
定義から`rfl`で一致する)**表現で書かれた新しいゴールに`rw`で適用
しようとすると、
```
Tactic `rewrite` failed: Did not find an occurrence of the pattern
  ...(psiGenStepResult K...m (psiGenSeq K...m)).pt...
in the target expression
  ...(psiGenSeq K...(m+1)).pt...
```
で失敗する。`rw`(`kabstract`)は**構文的な一致**しか見ないため、両者が
`rfl`で等しいという事実(`psiGenSeq`の定義展開だけで従う、`#34`で
既に確認済みの事実)を素通りしてしまう。

**直し方**: その箇所だけ`rw`ではなく、ゴールの形に**すでに一致する
型**を持つ項として引用する——`exact 定理名 引数...`、または(ゴールが
`f (…左辺…) = f (…右辺…)`のように**片側の関数適用として**書ける
場合は`congrArg f (定理名 引数...)`。どちらも`isDefEq`(kernelの
defeqチェック)で照合するため、構文的な一致は不要——`psiGenStepResult`
経由の型と`psiGenSeq`経由の型が(`rfl`で)同じである限り、そのまま
通る。今回は後者(`congrArg (principalUnitsQuotientEquiv K hπmax
(m+1) _)`で`reciprocityMapLimitCompat`を橋渡し)で解決した。

**教訓の一般化**: `#34`が教えた「個々の射影は`rfl`で軽いが、複数の
射影を組み合わせた**独立な型注釈**との照合は別問題として重い」と
表裏一体——今回は逆に、**個々の射影が`rfl`で一致する**という事実を
`rw`(構文照合)ではなく`exact`/`congrArg`(defeq照合)に**乗せ換える
だけ**で、新しい罠を踏むことなく橋渡しできた。`rw`で詰まったら、
まず「両辺は本当に`rfl`で一致するはずでは?」と疑い、`exact`/
`congrArg`に持ち替えてみるのが低コストな次の一手になる。

実例: `lean/ABC3/Found/PGC/LubinTateReciprocityMapLimit.lean`の
`reciprocityMapLimitFamily_step`(`m+1`の場合)。

**追記(2026-09-05、CorrHyp側の別の実例)**: 同じ病の別の症状として、
`rw [lemma1, lemma2]`のように**同じ補題を2回連鎖させる**場面でも
起きる——`lemma1`が`Y ⁻¹ᵁ (Z.basicOpen r)`の形を`Y`(定義上の展開形、
例えば`pullback X.hom toBaseK`)基準の項へ書き換えてしまい、2回目の
`rw`が要求する`α : C ⟶ (ExtF.obj X).left`(構文上の別名)との一致
チェックに失敗して「motive is not type correct」になる(`Scheme.
preimage_basicOpen`を`pullback.fst`・`α`の2段に適用する場面、
`ExtLimit.lean`の`piece_basicOpen_mul_eq`)。直し方は同じ:中間結果を
**明示的に型注釈した`have`**として確定させ(`(... : (ExtF.obj X).left.
Opens)`のように書きたい形を先に固定する)、`rw`ではなく`exact`(defeq
判定)で個別に閉じる——`rw`の自動連鎖を諦めて1段ずつ`have`で刻むのが
安全。

## 38. `Option.elim i A B`のような**依存する`match`族**に対する
`instance`は、`inferInstance`だけでは各枝で閉じない

`Option.elim i QC Q`(`i:none`なら`QC`、`i:some j`なら`Q j`)という
族に対して`∀ i, AddCommGroup (Option.elim i QC Q)`という`instance`を
`match`で

```lean
instance : ∀ i : Option α, AddCommGroup (Option.elim i QC Q)
  | none => inferInstance
  | some j => inferInstance
```

のように書くと、`none`の枝で `failed to synthesize instance of type
class AddCommGroup (none.elim QC Q)` のように**簡約前の形のまま**
探索に失敗する——`match`の各枝で`Option.elim`が定義通り`QC`/`Q j`に
簡約されることを、`inferInstance`(型クラス探索)は自動では認識
しない(探索は構文的な頭部記号でインスタンスを絞り込むため、
`Option.elim none QC Q`という頭部が`Option.elim`のままだと`QC`用の
インスタンスが候補に挙がらない)。直し方は、`show`で先に**簡約後の
型**へ変換してから`infer_instance`する:

```lean
instance : ∀ i : Option α, AddCommGroup (Option.elim i QC Q)
  | none => by show AddCommGroup QC; infer_instance
  | some j => by show AddCommGroup (Q j); infer_instance
```

`show`は`isDefEq`(kernelのdefeqチェック)で照合するので、
`Option.elim none QC Q`と`QC`が定義上等しいことは問題無く通る——
`#1`「instances透明度で型が合わない」・`#33`「型クラス**探索**にも
起きる」と同根だが、今回は「探索対象の型そのものが未簡約」という、
また別の症状。

実例: `lean/ABC3/Found/Falt1/KaehlerAux.lean`の
`falt1_optionElim_addCommGroup`/`falt1_optionElim_module`
(`pushoutKaehlerSplitStepOption`の`Option ι`出力と、成分ごとの
全射`φC`/`φ i`の像`QC`/`Q i`を束ねる場面)。

## 39. instance diamondの回避策として、`have`/`set`で独立に型注釈
するのではなく**呼び出し引数の位置に無名関数を直接書く**

`#1`「instances透明度で型が合わない」の一種だが、直し方が別角度:
`letI`で望む`Algebra`instanceを明示登録しても、その**後**に
`have φ : (…その instance が要求される型…) := ...` のように
`φ`を**独立に型注釈して`have`/`set`する**と、その型注釈自体の
elaborationが型クラス探索を独自に再実行し、`letI`で登録した
instanceではなく別の(グローバルな)instanceを見つけてしまうことが
ある——`letI`の登録は「後続のtermのinstance探索で優先される」とは
限らず、**新しい型注釈のスコープでの独立した探索**では負けうる。

直し方: `φ`を独立に`have`/`set`せず、**それを引数として渡す関数
呼び出しの引数位置に無名関数のまま直接書く**:

```lean
-- 悪い例(instance diamondで失敗する):
have φ : (期待される具体的な型) := fun i => 0
exact someTheorem ... φ ...

-- 良い例(呼び出し先の期待する型から直接推論させる):
exact someTheorem ... (fun i => 0) ...
```

こうすると`φ`の型は独立に決まらず、`someTheorem`の**仮引数の型**
(=呼び出し先が実際に要求する、正しいinstanceを含む型)から直接
推論されるため、独立した型注釈によるinstance探索の分岐がそもそも
発生しない。`#1`の対処法(`letI`で明示登録)がうまくいかない場面で、
まず試す価値がある軽量な代替。

実例: `lean/ABC3/Found/Falt1/KaehlerAux.lean`の
`falt1_pushoutKaehlerSplitStepOption_adjoinRoot_surjective_example`
(`RAlgOver.lift`由来の`Algebra`instanceと`AdjoinRoot.instAlgebra`の
diamond)。

## 40. `MvPolynomial ι (テンソル積の型) ⧸ I`の`HasQuotient`自動探索が
失敗する——`letI`で確定させるべきは**係数環自身**であって`MvPolynomial`
本体ではない(2026-09-05)

**症状**: `I : Ideal (MvPolynomial ι B)`(`B := A ⊗[ℚ] R.1`のような
テンソル積型)を書いた後、`MvPolynomial ι B ⧸ I`という式(`⧸`記法単体、
`descendPieceR`のような既存コードにも現れる形)が、
```
failed to synthesize instance of type class
  HasQuotient (MvPolynomial ι B) (Ideal (MvPolynomial ι B))
```
で失敗することがある——単純な`B := ℚ ⊗[ℚ] ℝ`のような最小例でも再現する
(`CommRing`型クラス自体ではなく`HasQuotient`のnotation展開が引く`Ideal`
の`Semiring`インスタンスの解決だけが、`MvPolynomial`本来の`CommRing`
インスタンスと非`defeq`に見える別の道(`AddMonoidAlgebra.semiring`
経由)を取ってしまうため)。`CommRing (MvPolynomial ι B)`という直接の
ゴールは`infer_instance`で普通に通る——`⧸`記法特有の失敗形。

**直し方**: `letI hCR : CommRing (MvPolynomial ι B) := inferInstance`
のように**`MvPolynomial`自身**のインスタンスを先に確定させても効果が
無い——効くのは**係数環`B`自身**のインスタンスを先に確定させること:
```lean
letI hCR : CommRing (A ⊗[ℚ] R.1) := inferInstance  -- ← Bの方
-- この後なら MvPolynomial ι (A ⊗[ℚ] R.1) ⧸ I が普通に通る
```
`descendPieceR`(`ExtLimit.lean`)がこのパターンを最初から採用していた
理由がこれで判明した——`#1`「instances透明度」系の変種だが、「どちらの
型のインスタンスを`letI`で確定させるべきか」を間違えると効かない、
という一段具体的な教訓。

実例: `lean/ABC3/Found/CorrHyp/FieldLimit.lean`の
`exists_fg_subalgebra_tensor_quotientMvPolynomial_lift`。

## 41. 型注釈中の**無名の`let`(インスタンス以外)**は、証明本体で
`intro`しないと後続の`∀`/`∃`束縛と名前が衝突する——`letI`/`haveI`
(インスタンス用)とは違って**自動ではゼータ簡約されない**(2026-09-05)

**症状**: `letI := ...`(インスタンス)は型注釈の中に置いておけば証明
本体で`intro`しなくても(自動的にゼータ簡約されて)そのまま使える
(`descendPieceR_toBase`等、既存コードの標準パターン)。ところが**同じ
型注釈の中に、インスタンスではない普通の値を束縛する無名の`let`**
(例: `let n := Algebra.Presentation.ofFinitePresentationVars ...`)を
置き、その**後**に`∀ (e : ... n ...), ...`のように`n`を参照する束縛を
続けると、証明側でこの`let`は自動簡約**されない**——`intro e`を最初に
呼ぶと、実際には`∀e`ではなく**この`let n`を`e`という名前で消費して
しまう**(束縛の順序どおりに`intro`は`let`も1個ずつ数える)。結果、
`e`が本来の型ではなく`let`の値の型(この例では`ℕ`)を持つことになり、
`e.symm`のような呼び出しが`Nat.symm`を探して失敗するという、原因が
非常に分かりにくいエラーになる(`intro`自体はエラーを出さず`成功`
してしまうため、型不一致のエラーが**ずっと後**の行に出る)。

**直し方**: 型注釈中に置いた無名の`let`は、宣言された**順番どおりに
全部**`intro`で明示的に消費する——`let n := ...`・`let I := ...`の後に
`∀ e : ...`があるなら、`intro n I e`と1回で3つとも(または`intro n;
intro I; intro e`と分けて)消費する。証明側で`set R := ...`のような
別名を与えていても、それは**別の`have`/`letI`**であって型注釈中の
`let`の`intro`を代替しない。

**波及**: 「巨大な依存型を`∀`の中で複数回書くと`maxHeartbeats`が
(4000000でも)尽きる」という、以前(`#40`寄りの文脈)観測した現象は、
実際にはこの`intro`忘れが真因だった可能性が高いと判明した——`let`で
式を共有した上で`intro`を正しい個数・順序で行えば、`maxHeartbeats
4000000`(実測82秒)で普通に通った。

実例: `lean/ABC3/Found/CorrHyp/ExtLimit.lean`の
`exists_piece_basicOpen_R_lift`。

## 42. `open ... in`・`set_option ... in`は`/-- docstring -/`より
**前**に置く——後に置くと構文エラーになる(2026-09-05)

**症状**: `/-- 説明文 -/`の直後に`open scoped X in`を置き、その次の行に
`theorem ...`を書くと、
```
error: unexpected token 'open'; expected 'lemma'
```
という構文エラーになる——`open`という単語自体は正しいのに、パーサーが
`lemma`(または`theorem`/`def`)を期待している場所に来てしまう。

**直し方**: 修飾子コマンド(`open ... in`・`set_option ... in`)は
**docstringより前**に置く:
```lean
-- 良い例
open scoped TensorProduct in
/-- 説明文 -/
theorem foo ... := ...

-- 悪い例(構文エラー)
/-- 説明文 -/
open scoped TensorProduct in
theorem foo ... := ...
```
`docstring`は直後の宣言キーワード(`theorem`/`def`/`lemma`等)に直接
くっつく必要があり、間に`open ... in`のような修飾子コマンドを挟めない
——ファイル内の既存コード(`piece_le_of_le`等)はすべて「修飾子→
docstring→宣言」の順で書かれており、今回はその逆順にしてしまった。

実例: `lean/ABC3/Found/CorrHyp/ExtLimit.lean`の
`descendPieceR_localization_isOpenImmersion`。

## 43. `IsLocalization.ringEquivOfRingEquiv`等を具体的な型(テンソル積
など)のまま直接呼ぶとinstance diamondに当たる——**先に抽象的な型変数
の補題として切り出してから代入する**と回避できる(2026-09-05)

**症状**: `e2 : (A⊗[R]B) ≃+* C`(具体的なテンソル積の型)を持っていて
`IsLocalization.ringEquivOfRingEquiv S Q e2 proof`(局所化同士の同型を
作る、`Localization.Away`に使う定番)を直接呼ぶと、
```
Application type mismatch: The argument e2 has type
  RingEquiv ... Algebra.TensorProduct.instMul instDistribOfSemiring.toMul ...
but is expected to have type
  RingEquiv ... instDistribOfSemiring.toMul instDistribOfSemiring.toMul ...
```
という、`A⊗[R]B`の`Mul`/`Add`インスタンスが(`Algebra.TensorProduct.
instMul`経由か`instDistribOfSemiring.toMul`経由かで)非`defeq`に見える
instance diamondになる——最小の反例(`A⊗[R]B`を含む式で`Localization.
Away`を直接構成しようとする式)でも再現する、`#1`「instances透明度」
系の中でも特に頑固な変種。`letI`での事前登録(`#1`の定石)を`Mul`側・
`Add`側それぞれに何度試しても解消しない。

**直し方**: `IsLocalization.ringEquivOfRingEquiv`(や同種のAPI)を、
**具体的な型(テンソル積等)のまま直接使わず**、まず**抽象的な
`CommRing`型変数**だけを引数に取る小さな補題として**別立てで**用意
してから、その補題へ具体的な型を**代入して使う**:
```lean
-- 抽象化した補題(具体的な型を一切知らない、diamondが起きない)
theorem ringEquiv_localization_of_apply_eq (A B : Type) [CommRing A] [CommRing B]
    (e : A ≃+* B) (a : A) (b : B) (hab : e a = b) :
    Nonempty (Localization.Away a ≃+* Localization.Away b) :=
  ⟨IsLocalization.ringEquivOfRingEquiv (Localization.Away a) (Localization.Away b) e
    (hab ▸ Submonoid.map_powers e.toMonoidHom a)⟩

-- 使う側: A・Bに具体的な型(テンソル積など)を代入するだけ、diamond無し
obtain ⟨e3⟩ := ringEquiv_localization_of_apply_eq (A⊗[R]B) C e2 a b hab
```
`letI`で「同じ場所に」インスタンスを事前登録するのとは**違う対処**
——`letI`は型を固定したまま探索順序だけを変えようとするが、抽象化は
**型変数を経由させることで、具体的な型の内部構造(`Algebra.TensorProduct.
instMul`のような複合インスタンス)自体をelaboratorに一切見せない**、
という一段違うレベルの回避策。`Away`局所化に限らず、`IsLocalization`
系のAPI全般で同じ手が効く可能性がある。

実例: `lean/ABC3/Found/CorrHyp/FieldLimit.lean`の
`ringEquiv_localization_of_apply_eq`・`exists_ringEquiv_localization_of_eq`。

## 44. `Algebra.Etale`/`Module.Free`/`Module.Finite`のAlgEquiv移送で
instance diamondに落ちる時は、**`RingHom.Etale`(bare ring homの性質)
のレベルまで降りて`IsLocalization.ringHom_ext`で押し切る**と迂回できる
(2026-09-05)

**症状**: `awayAlgebra`のような`(f).toAlgebra`で人工的に配線した
`Algebra`インスタンスに対して、自然に存在するはずの`AlgEquiv`
(例: `R ≃ₐ[R] Localization.Away 1`・`Fin2→R ≃ₐ[Fin2→R] Localization.Away
(algebraMap R (Fin2→R) 1)`)を使って`Algebra.Etale`を移送しようとすると、
`▸`・`convert`・`AlgEquiv`合成のどこかで人工instanceと標準instanceが
非`defeq`に見えて詰まる(`#1`・`#43`と同系統だが`Algebra.Etale`は
`Prop`なので`convert`の残り目標が余計に混乱しやすい)。

**直し方**: `Algebra`インスタンス・`AlgEquiv`を経由せず、**bare
`RingHom`の性質`RingHom.Etale`(`Mathlib.RingTheory.RingHom.Etale`)の
レベルで組み立てる**——`RingHom.Etale`はinstanceを一切参照しない
「`f : R →+* S`の性質」なので diamond がそもそも起こらない:
```lean
-- 1. 目的の環準同型(`awayAlgebra`が使う`Localization.awayMap ...`)を、
--    「全単射(=同型) ∘ 素の algebraMap ∘ 全単射(=同型)」に分解する式を、
--    AlgEquiv ではなく RingEquiv.ofBijective + RingHom の等式として、
--    IsLocalization.ringHom_ext(局所化からの環準同型の一意性)で示す:
have heqmap : Localization.awayMap (algebraMap R (Fin2→R)) 1
    = ιB.toRingHom.comp ((algebraMap R (Fin2→R)).comp ιR.symm.toRingHom) := by
  apply IsLocalization.ringHom_ext (Submonoid.powers (1:R))
  ext x; ...  -- unfold Localization.awayMap IsLocalization.Away.map; rw [IsLocalization.map_eq]
-- 2. 各ピースをRingHom.Etaleで示し、stableUnderComposition/of_bijectiveで貼る:
have hgEtale : RingHom.Etale (algebraMap R (Fin2→R)) := RingHom.etale_algebraMap.mpr etale_fin2
have hιEtale : RingHom.Etale ιR.symm.toRingHom := RingHom.Etale.of_bijective ιR.symm.bijective
have hcomp : RingHom.Etale (... .comp ...) :=
  RingHom.Etale.stableUnderComposition ιR.symm.toRingHom (algebraMap R (Fin2→R)) hιEtale hgEtale
-- 3. 最後に heqmap で書き換えてから Algebra インスタンスへ戻す:
rw [← heqmap] at hfullEtale
exact RingHom.Etale.toAlgebra hfullEtale
```
`Module.Free`/`Module.Finite`側は`AlgEquiv`ではなく**半線形同値**
(`≃ₛₗ[σ]`、`σ`は`ιR.toRingHom`のような具体的なRingHom)を`heqmap`から
組み立て、`Module.Free.of_equiv`(半線形移送)・`Module.Finite.
of_equiv_equiv`(`≃ₗ`を経由しない、2つのRingEquivと可換四角形からの
移送——`Module.Finite.equiv`は同じ環上の`≃ₗ`しか受け付けないため使えない)
に渡す。

**注意点2つ**: (a) `(e : R →+* S) = algebraMap R S := by ext x; exact
e.commutes x`はドメインが`Fin2→R`のような`Pi`型だと`ext`が`Pi.single`
方向に過分解して失敗する——`RingHom.ext (fun x => e.commutes x)`を使う。
(b) `letI`が**定理の型**にしかない場合、証明本体の中で`Module.Free.
of_equiv`等がinstance探索に失敗する——証明の最初の行で`letI`を
再宣言する(`awayScalarTower`の証明が既にこの形)。

実例: `lean/ABC3/Found/Falt1/AlmostEtale.lean`の
`awayOne_fin2_etale`・`awayOne_fin2_freeFinite`(`Definition 2.1`の
witness、`A:=R`・`B:=Fin2→R`・`p:=1`の場合の条件(i)の`Etale`・
`Free`・`Finite`部分)。

## 45. `Exists.choose`で非構成的に定義された的の元(`elem`等)への到達を
示す時、**その元を計算しようとせず、経由する写像の全射性だけを示す**
と、`Exists.choose`の中身を一切知らずに済む(2026-09-05)

**症状**: `Algebra.FormallyUnramified.elem R S : S⊗[R]S`(diagonal
idempotent)のような、`(iff_exists_tensorProduct.mp inferInstance)
.choose`で定義された元は、**mathlibに具体的な計算式が無い**
(存在と一意特徴付けの性質——`lmul_elem`・`one_tmul_sub_tmul_one_mul_
elem`——しか使えない)。`∃e, f e = elem R' S'`(`f`は何らかの環準同型
から誘導されるテンソル積上の写像)を示したい時、素朴には「`elem R S`
を明示的な元(`e1⊗e1+e2⊗e2`等)だと特定してから`f`で送る」という
経路を考えがちだが、これは`elem R' S'`が「たまたまその元と等しい」
ことを`Exists.choose`の非構成性のせいで**別途、一意性補題を経由して
証明する**必要が生じ、大掛かりになる(一意性補題自体もmathlibに無い)。

**直し方**: 目的の元(`elem R' S'`)を計算・特定しようとせず、
**`f`が全射であることだけを示す**——全射なら`elem R' S'`が何であれ
(値を一切知らなくても)`∃e, f e = elem R' S'`は`f`の全射性の定義
そのものから直ちに従う。全射性は`TensorProduct.induction_on`
(`zero`・`tmul`・`add`の3ケース)で、`tmul`ケースだけ「土台の環準同型
`f0:B→Bp`が全射」という**遥かに単純な事実**(今回は`p`が単元である
ことからの`IsLocalization.atUnit`の全単射性)に帰着させれば良い。
```lean
-- f0 が全射なら、f0を両成分に施すだけの写像も全射
theorem diagonalCompare_surjective_of_algebraMap_surjective ... (hsurj : Function.Surjective f0) :
    Function.Surjective (diagonalCompare p) := by
  intro z
  refine TensorProduct.induction_on z ⟨0, map_zero _⟩ ?_ ?_
  · intro x y
    obtain ⟨b1, hb1⟩ := hsurj x; obtain ⟨b2, hb2⟩ := hsurj y
    exact ⟨b1 ⊗ₜ b2, by rw [diagonalCompare_tmul, hb1, hb2]⟩
  · rintro x y ⟨ex, hex⟩ ⟨ey, hey⟩
    exact ⟨ex + ey, by rw [map_add, hex, hey]⟩
-- 使う側: elem の値を一切知らずに existence が出る
obtain ⟨e, he⟩ := hsurj (Algebra.FormallyUnramified.elem Ap Bp)
```
注意点: `letI`/`haveI`を型シグネチャに埋め込んだ定理(`diagonalCompare`
本体等)を`Function.Surjective (diagonalCompare p)`のように**引数無し
で関数として渡す**と、`diagonalCompare`の暗黙引数`{A B}`のうち`p`の
型からだけでは決まらない方(`B`)が推論できず`typeclass instance
problem is stuck`になる——`diagonalCompare (A := A) (B := B) p`と
**明示的に埋める**必要がある(`#1`の変種、関数を値として渡す時は
暗黙引数を推論に任せきらない)。

もう1点: `Algebra.FormallyUnramified.elem Ap Bp`のような`[Formally
Unramified][EssFiniteType]`要求の項を**定理の型シグネチャ内**で使う
時、`Algebra.Etale`・`Module.Finite`からこれらを`infer_instance`で
導出できても、その導出元の`haveI`は**シグネチャ自身の中で**(証明
本体の中ではなく)`letI := ...; haveI := ...; <本体の型>`の形で先に
宣言しておく必要がある——証明本体内の`haveI`は「型を検査する時点」
より後に実行されるため、型そのものの中の`elem`の要求するインスタンス
解決には間に合わない(`IsAlmostEtaleCovering`・`diagonalCompare`
自身が既にこの`letI`/`haveI`前置パターンを使っている理由)。

実例: `lean/ABC3/Found/Falt1/AlmostEtale.lean`の
`diagonalCompare_tmul`・`diagonalCompare_surjective_of_algebraMap_
surjective`・`awayOne_fin2_idempotent`・`awayOne_fin2_
isAlmostEtaleCovering`(`Definition 2.1`の non-vacuous witness、
条件(iii)・最終組み立て)。

## 46. 2つの `{A B : Type*}` が「同じ universe」だと思い込むと、
`RingHom.StableUnderComposition` 系の補題(`{R S T : Type u}` と**単一の**
universe 変数を要求)が意味不明な type mismatch で失敗する
(2026-09-05)

**症状**: `#44`の戦略(`RingHom.Etale.stableUnderComposition f g hf hg`)
を、具体的な型(`R`・`Fin 2 → R`)ではなく**抽象的な2つの型変数**
`{A B : Type*}`に一般化しようとすると、`f : ιA.symm.toRingHom`・
`g := algebraMap A B`を渡しただけの、具体版と全く同じ形の呼び出しが
```
Type mismatch
  RingHom.Etale.stableUnderComposition ιA.symm.toRingHom ?m.592 hιAinvEtale ?m.593
has type
  (fun {R S} [CommRing R] [CommRing S] => RingHom.Etale) (RingHom.comp ?m.592 ιA.symm.toRingHom)
but is expected to have type
  ((algebraMap A B).comp ιA.symm.toRingHom).Etale
```
という、`g`(`?m.592`)が最後まで決まらない type mismatch で失敗する
——`(f := ...)`/`(g := ...)`という named argument をつけても変わらない。

**原因**: `RingHom.StableUnderComposition`(`RingHom.Etale.
stableUnderComposition`の型)は`∀ ⦃R S T : Type u⦄ ...`と**単一の
universe 変数 `u`** を要求する。`{A B : Type*}`と素朴に書くと、Lean は
`A`・`B`にそれぞれ**別々の** universe metavariable(`Type u_1`・
`Type u_2`)を割り当てる——`stableUnderComposition`を`Localization.
Away 1 : Type u_1`(`A`由来)と`B : Type u_2`の間の合成に使おうとした
瞬間、`u_1 ≠ u_2`の可能性を排除できず、統一に失敗する(具体的な型
`R`・`Fin 2 → R`ではこの2つが**たまたま同じ**universeだったため、
一般化するまで症状が出なかった)。

**直し方**: `{A B : Type*}`をやめ、`universe u`を宣言してから
`{A B : Type u}`と**明示的に同じ**universe変数を使う。
```lean
universe u
theorem foo {A B : Type u} [CommRing A] [CommRing B] [Algebra A B] ... := by
  ...
  exact RingHom.Etale.stableUnderComposition f g hf hg  -- これで通る
```
`AlgHom`/`RingHom`の合成安定性(`StableUnderComposition`)系の補題
全般に共通する注意点——`Algebra.TensorProduct.lift`等、同一universeを
要求しない多くのAPIでは`Type*`のままで問題ない。

実例: `lean/ABC3/Found/Falt1/AlmostEtale.lean`の`awayOne_etale_of_etale`
(`Definition 2.1`のwitnessを`B:=Fin2→R`から「`Etale`・`Finite`・`Free`
な任意の`B`」へ一般化、`B:=A`(恒等拡大)の場合を含む)。

## 48. `Exists.choose`で定義された元(`elem`等)を「値を計算せず一意性
だけで」他の環へ自然に移す時は、**部分環の生成元閉包**(`+`・`*`・
逆元)で議論を単項生成から全体へ拡張する(2026-09-05)

**症状**: `#45`の「全射性だけで押し切る」トリックは `p` が単元の
時だけ効く(`algebraMap B Bp` が全単射になるのは局所化が退化する
`p`単元の場合のみ)。`p` が真の非単元素元の一般の場合、`elem A B` の
局所化での像(`diagonalCompare p (elem A B)`)が `elem Ap Bp` に
一致することを示すには、**`elem` の定義性質(1)を`Bp`の全ての元
`s'`について**確認する必要があるが、`s'`は`f0(B)`(`f0:B→Bp`、一般に
全射でない)の像とは限らない。

**直し方**: `Z:={s'∈Bp | (1⊗s'-s'⊗1)*t=0}`という述語を考え、これが
**部分環**であること(`0`・`+`・`*`で閉じる、直接の代数計算——
`(a⊗1)*[(1⊗b)-(b⊗1)]+[(1⊗a)-(a⊗1)]*(b⊗1) = 1⊗(ab)-(ab)⊗1`という
恒等式が鍵)、かつ **`Away`局所化の生成元 `π:=f0(p)` の逆元も`Z`に
入る**こと(`π`自身は`elem`の定義性質(1)から`Z`に入るとわかっている
——`(π⊗π)`を掛けてから可逆性でキャンセルする論法、`π`が単元なのは
`Away`局所化の定義そのものから)を示す。最後に`IsLocalization.
Away.surj`(`Bp`の任意の元は`f0(a)*π⁻ⁿ`の形)で「`f0(B)`と`π⁻¹`から
`Z`が生成される」ことと「`Z`が部分環」を組み合わせて、`Z=Bp`(全射性
ではなく**全域性**)を結論する。「一意性補題」(`t,t'`が両方定義性質を
満たすなら`t=t'`——`t*(1⊗s-s⊗1)`の吸収性質`x*t=(μx⊗1)*t`から
`(1-t')*t=0`⟹`t=t'*t`、対称に`t'=t*t'`、可換性で`t=t'`)と組み合わせる
ことで、**`elem`の値を一切計算せずに**局所化での自然性
(`diagonalCompare p (elem A B) = elem Ap Bp`、`p`が単元かどうかに
一切依存しない)を証明できる。

実例: `lean/ABC3/Found/Falt1/AlmostEtale.lean`の`elem_unique_of_props`
(一意性)・`Zclosed_add`/`Zclosed_mul`/`Zclosed_inv`(部分環の生成元
閉包)・`diagonalCompare_elem_eq`(自然性、任意の`p`)・
`isAlmostEtaleCovering_of_etale_general`(最終組み立て、`Definition
2.1`のwitnessを`p`が単元でない**真に非退化な**場合まで一般化)。

## 47. `maxHeartbeats`の`whnf`タイムアウトは「構造的に不可能」を意味
しない——桁を上げて気長に待つだけで通ることがある(2026-09-05)

**症状**: `pieceAlgebra`等`CorrHyp`固有の巨大な足場(`letI`が10個
以上)を伴う証明が、`set_option maxHeartbeats 4000000`(既存コードの
標準的な上限)でも`(deterministic) timeout at whnf`になる。この症状
だけを見ると「この組み合わせ自体が構造的に無理(instance diamond等)」
と判断してしまいがちだが、**実際には単に計算資源が足りていないだけ**
のことがある。

**確認の仕方**: `maxHeartbeats`を大幅に(例えば10倍、`40000000`)
上げて、`mcp__abc3-lean__lean_check`をバックグラウンドで実行し
(`timeoutSeconds`を大きく、または既定の120秒超過で自動的にバック
グラウンド化されるのに任せ)、気長に待つ。今回の実例では`4000000`
で`whnf`タイムアウトしていた配線が、`40000000`で**227秒**かけて
無事通った——「証明できない」のではなく「もっと時間がかかる」だけ
だった。

**教訓**: `letI`の数が多い巨大な証明でタイムアウトに当たったとき、
最初に疑うべきは(1)`intro`忘れ(`#41`)、(2)欠けているインスタンス
(`letI`の追加漏れ)、(3)`algebraMap`記法の曖昧さ(`#1`系)——これらを
1つずつ潰してもなお`whnf`タイムアウトが残るなら、**設計を諦める前に
まず`maxHeartbeats`を10倍にして試す**のが低コストな次の一手になる。
ただし個々の`lean_check`呼び出しが数分かかるようになるため、最終的に
`.lean`ファイルへ書く際は`lake build`(こちらは`mathlib`用のグローバル
設定に従うため、ファイル冒頭の`set_option maxHeartbeats N in`を宣言に
直接付ける必要がある)でも同じ上限を明示することを忘れないこと。

実例: `lean/ABC3/Found/CorrHyp/ExtLimit.lean`の
`exists_descendPieceR_localization_baseChange`(`maxHeartbeats
40000000`、227秒)。

## 49. `IsLocalization`インスタンスを環同型に沿って移送したいとき、
「両立する`Algebra`構造を`e∘algebraMap`として定義する」だけで
`AlgEquiv`化でき、そのまま`IsLocalization.isLocalization_of_algEquiv`
へ渡せる(2026-09-05)

**状況**: `IsLocalization M S`が分かっていて、別の環`P`と環同型
`e : S ≃+* P`があるとき、「`P`もまた`M`による局所化である」ことを
示したい。`P`にはまだ`R`上の`Algebra`構造が無い(または、あっても
`e`と両立する保証が無い)ことが多い。

**やってはいけないこと**: `e`と既存の`Algebra R P`インスタンス(もし
既にあれば)が両立することを別途証明しようとする、あるいは
`IsLocalization`の定義(3条件:単元性・全射性・核の同値関係)を素朴に
展開して`e`越しに1つずつ確認する——長く、事故りやすい。

**正しいやり方**: `letI : Algebra R P := (e.toRingHom.comp
(algebraMap R S)).toAlgebra`で`P`の`Algebra R P`構造を**`e`を通した
ものとして定義**すると、`algebraMap R P`の定義がまさに`e∘algebraMap
R S`になるので、両立条件`∀r, e(algebraMap R S r) = algebraMap R P r`
は`rfl`で終わる。あとは`AlgEquiv.ofRingEquiv (f:=e) he`(mathlib、
両立条件から`RingEquiv`を`AlgEquiv`へ格上げ)+`IsLocalization.
isLocalization_of_algEquiv`(mathlib、`AlgEquiv`越しに`IsLocalization`
インスタンスを移送)を合成するだけ。

```lean
theorem isLocalization_of_ringEquiv_transport (R S P : Type) [CommRing R] [CommRing S] [CommRing P]
    (M : Submonoid R) [Algebra R S] [IsLocalization M S] (e : S ≃+* P) :
    letI : Algebra R P := (e.toRingHom.comp (algebraMap R S)).toAlgebra
    IsLocalization M P := by
  letI : Algebra R P := (e.toRingHom.comp (algebraMap R S)).toAlgebra
  have he : ∀ r, e (algebraMap R S r) = algebraMap R P r := fun r => rfl
  exact IsLocalization.isLocalization_of_algEquiv M (AlgEquiv.ofRingEquiv (f := e) he)
```

**教訓**: 「環同型越しにインスタンスを移送したい」系の問題は、まず
「移送先の`Algebra`/加群構造を環同型の**定義そのものとして**構成
できないか」を考えると、両立条件が`rfl`で落ちて一気に軽くなる。
`decl-index.txt`の検索結果(`isLocalization_iff_of_ringEquiv`等)が
実際には見つからない(namespace違い、または版違い)ことがあるので、
見つからないときは`exact?`で型から直接引く(`IsLocalization.
isLocalization_of_algEquiv`はこの方法で発見した)。

実例: `lean/ABC3/Found/CorrHyp/FieldLimit.lean`の
`isLocalization_of_ringEquiv_transport`。

## 50. `B⊗_AB`から`Bp⊗_ApBp`への比較写像の単射性・「後続の仮説の型が
参照している仮説」を`obtain`/`unfold ... at`で壊す罠・`AlgHom`越しの
`map_smul`が`map_smul_of_tower`を要る場面、の3点セット(2026-09-05)

**(a) 局所化を2回経由したテンソル積を1回に戻す同型**: `M₁,M₂`が
`R`-代数`A`(`IsLocalization S A`)上の加群のとき、`TensorProduct A M₁
M₂ ≃ₗ[A] TensorProduct R M₁ M₂`という「`A`上のテンソルは、実は`R`上の
テンソルと(`M₁,M₂`が既に局所化されている分)同じ」という同型が
`IsLocalization.moduleTensorEquiv (S) (A) (M₁) (M₂)`として存在する
(pure tensorを pure tensor に送るだけ、`exact?`で`e(x⊗y)=x⊗y`が
確認できる)。`TensorProduct.AlgebraTensorModule.cancelBaseChange`や
`IsLocalization.Away.tensorRightEquiv`/`tensorEquiv`ではこの形の
同型は得られない(仮定の噛み合わせ・向きが合わない)——この形が
必要なら`moduleTensorEquiv`を先に探す。

これを使うと、環準同型`f : B⊗_AB→Bp⊗_ApBp`(`Bp,Ap`は`p`による
局所化)が`φ:=(rTensor Bp f₀)∘(lTensor B f₀)`(`f₀:B→Bp`をテンソル積の
両成分に適用するだけ)と`e:=moduleTensorEquiv`の合成`e.symm∘φ`に
分解できるとき、`φ`の単射性(`Module.Flat`の`rTensor`/`lTensor`
`preserves_injective_linearMap`を2回、`f₀`の単射性から)だけで`f`の
単射性が言える——`Module.Flat A B`は`Module.Free A B`から`inferInstance`、
`Module.Flat A Bp`は`IsLocalization.flat`+`Module.Flat.trans`。

**(b) `obtain ⟨...⟩ := h`(または`unfold ... at h`)は、後続の仮説の
型が`h`を参照していると、その仮説を静かに消す**: `theorem foo (h :
∃ x, P x) (hw : (h.2.1 を使う型)) : ... := by obtain ⟨x, hx⟩ := h; ...`
のような形で、証明中に`h`を`obtain`で分解すると、`hw`の型が(古い)`h`
という自由変数を参照しているため、Lean は`hw`を revert してから`h`を
分解する。この revert-and-reintroduce が期待通りに`hw`という名前の
まま戻ってくると思い込むと、`lake build`で`Unknown identifier hw`や
`unsolved goals`(その少し手前の行)として現れる——**しかも`lean_check`
(REPL)では検出できず`lake build`で初めて再現することがある**(今回
実際に踏んだ:REPLでは`OK`、ファイルへ書いて`lake build`したら失敗)。
**直し方**: `h`を`obtain`で壊さず、`haveI := h.1`・`haveI := h.2.1`
のように**射影(`.1`/`.2.1`等)をそのまま`haveI`/`have`に渡す**(`h`
自体は生かしたまま)。これなら`hw`の型が参照している`h`はそのまま
残るので何も壊れない。`h`の型が`def ... : Prop := ∃ ...`のような
略記(`IsAlmostEtaleCovering`等)でも、射影記法は定義展開して普通に
使える。**教訓: REPLでの`OK`は`lake build`の代わりにならない
(このファイルの「作業効率」条項の通り)——依存する仮説の型が絡む
局面では特に、最終確認を省略しない。**

**(c) `AlgHom R S`(`R`線形にしか型付けされていない)を`A`-scalar
(`IsScalarTower A R S`で`R`より"下"にある環)越しに`map_smul`したい
とき、生の`map_smul`は`MulActionHomClass`探索に失敗する**——`AlgHom`
は`R`-linearであることしか型に持たないため。`(f : M →ₐ[R] N).
toLinearMap.map_smul_of_tower (a : A) (x : M)`(`LinearMap.
map_smul_of_tower`、`IsScalarTower A R M`・`IsScalarTower A R N`を
使う版)に切り替えると通る。ただし`rw`で直接使おうとすると`(f :
M→ₐ[R] N) x`(AlgHomのcoe適用)と`f.toLinearMap x`(LinearMapのcoe
適用)が構文的に一致せずパターンが見つからないことがある——`have hpull
:= f.toLinearMap.map_smul_of_tower a x`を単独の項として取り出し、
`rw [show (f x) = a•(f y) from hpull]`のように`show ... from`で
包んで目的の構文形に揃えてから使うと安定する。

実例: `lean/ABC3/Found/Falt1/AlmostEtale.lean`の
`diagonalCompare_injective`(a)・`almost_swap_annihilate`(b)・
`almost_swap_augment`(c)。

## 51. `Ext X Y n`の`R`-加群構造は`Ext/Linear.lean`(別 import)にある。
かつ`@Foo.bar`が「unknown identifier」ではなく**インスタンス探索の
タイムアウト**を返すせいで、名前の存在確認自体が誤誘導される
(2026-09-05)

**状況**: `CategoryTheory.Abelian.Ext`(導来圏経由の一般 Ext)で
`r • e = 0`(`r`は係数環の元、`e : Ext X Y n`)のような**加群構造を
使う**主張を書こうとすると、`Module R (Ext X Y n)`のインスタンス
探索が失敗する。`Linear R C`(例:`ModuleCat.instLinear`で
`Linear T (ModuleCat T)`)も`HasExt C`も個別には`infer_instance`で
出るのに、`Module R (Ext X Y n)`だけが出ない。

**原因**: その`Module`インスタンスは
`Mathlib/Algebra/Homology/DerivedCategory/Ext/Linear.lean`にある
無名インスタンスで、`Ext/Basic.lean`や`Ext/EnoughProjectives.lean`
からは**推移的に import されない**。`import Mathlib.Algebra.Homology.
DerivedCategory.Ext.Linear`を1行足すだけで解決する。同ファイルには
`Ext.smul_comp`・`Ext.comp_smul`・`Ext.mk₀_smul`・
`Ext.smul_eq_comp_mk₀`・`Ext.linearEquiv₀`も入っている。

**なぜ気づきにくいか(ここが本題)**: import されていないことを
確認しようとして`#check @CategoryTheory.Abelian.Ext.mk₀_smul`と
書くと、**「unknown identifier」ではなく**
`Localization.HasSmallLocalizedHom`のインスタンス探索
タイムアウトが返る。`Ext`が`def`なので`Ext.mk₀_smul`が generalized
field notation(`Ext`を適用してからフィールドを取る)として解釈され、
`Ext`の暗黙引数の解決に入ってしまうため。実際、**存在しない名前**
`@CategoryTheory.Abelian.Ext.this_name_does_not_exist`でも
**同じエラー**が出る——つまりこのエラーは名前の有無について何も
語らない。

**確認の仕方**: 名前空間が`def`と衝突する場合は、`#check`ではなく
`have h := @Foo.bar`(タクティクブロック内)や、そもそも
`node tools/decl-index.mjs --mathlib`で作った`.cache/mathlib-index.txt`
を grep して**宣言がどのファイルにあるか**を見る(この索引は
ファイルパス付きなので、import の要否が同時に分かる)。今回は
`grep "Ext/Linear.lean" .cache/mathlib-index.txt`で一発だった。

**ついでに**: `set_option maxHeartbeats N` は**インスタンス探索の
上限を上げない**(そちらは`synthInstance.maxHeartbeats`、既定
20000)。`Ext`まわりは`HasSmallLocalizedHom`の探索が重く、
抽象的な環だと既定を超えることがあるので、
`set_option synthInstance.maxHeartbeats 400000 in`も併せて覚えておく
(ただし「探索が速く失敗する」場合は上限ではなくインスタンス不在
——今回のように import 漏れを疑う)。

実例: `lean/ABC3/Found/Falt1/AlmostEtale.lean`の
`ext_smul_eq_zero_of_almost_split`・`hochschild_ext_almost_zero`。

## 52. 商環(`... ⧸ I`)に`RingEquiv`越しの`Algebra`構造を載せる作戦は、
商環が**既に持っている**`Submodule.Quotient.instSMul'`に負けて破綻する
——最初から`AlgEquiv`(`≃ₐ[基底]`)を作るのが正解(2026-09-05)

**症状(3段のドミノ)**: `e : S ≃+* M`(`M`は`MvPolynomial ι B ⧸ J`)から
`letI : Algebra Q M := (e.toRingHom.comp (algebraMap Q S)).toAlgebra`と
インスタンスを移送し、さらに`letI : Algebra B M := ((algebraMap Q M).comp
(algebraMap B Q)).toAlgebra`と合成して`IsScalarTower B Q M`を`IsScalarTower.
of_algebraMap_eq (fun _ => rfl)`で出そうとすると、`lake build`のたびに
エラーが**1行ずつ先へ進む**:

1. `Algebra Q M`を作る行で`M`の`Semiring`が`Ideal.Quotient.semiring`と
   `CommSemiring.toSemiring`の2経路に割れる → `set M := ...`の**直後**に
   `letI hCRM : CommRing M := inferInstance`を置くと解消(`#40`と同型)。
2. 次に`Q`側で同じ割れ方をする → `set Q := ...`の直後にも
   `letI hCRQ : CommRing Q := inferInstance`が要る。
3. そして`IsScalarTower B Q M`の**型そのもの**が
   `@IsScalarTower B Q M (Submodule.Quotient.instSMul' I) algQM.toSMul
   (Submodule.Quotient.instSMul' J)`と表示され、`of_algebraMap_eq`が返す
   `Algebra.toSMul`3本組と合わない。

**根本原因**: 商環`MvPolynomial ι B ⧸ J`は、`B`が係数環である以上
**自前の**`SMul B _`(`Submodule.Quotient.instSMul'`)を持っており、
`SMul`のインスタンス探索ではこちらが勝つ。こちらは「係数への作用を
商へ降ろしたもの」であり、`e`越しに移送した`Algebra B M`(局所化の
構造を経由する別の写像)とは**構文的にも定義的にも一致しない**。
ローカルの`letI`で上書きしても`SMul`の層で負けるので勝てない。

**正しいやり方**: 移送で誤魔化さず、**最初から`AlgEquiv`を作る**。
つまり`S ≃+* M`ではなく`S ≃ₐ[B] M`を構成する。mathlibの部品は
たいてい`AlgEquiv`版が揃っている:
`IsLocalization.Away.mvPolynomialQuotientEquiv`(`≃ₐ[R]`)・
`MvPolynomial.quotientEquivQuotientMvPolynomial`(`≃ₐ[R]`)・
`DoubleQuot.quotQuotEquivQuotSupₐ`(`≃ₐ[R]`、`R`は明示引数)・
`Ideal.quotientEquivAlg`(`Ideal.quotientEquiv`の`≃ₐ`版)・
`MvPolynomial.sumAlgEquiv`(元から`≃ₐ`)。基底が途中で変わる箇所は
`AlgEquiv.restrictScalars`で下の基底へ落とす。

**教訓**: 商環・局所化のように「既定のインスタンスを自前で持っている」
型へインスタンスを後付け移送するのは、`Semiring`の割れ(1)(2)を潰しても
`SMul`の層(3)で必ず詰む。`≃+*`を作った時点で「あとで`Algebra`を移送
すればいい」と考えないこと——**基底を決めて`≃ₐ`で作る**のが唯一の
安定路線。

実例: `lean/ABC3/Found/CorrHyp/ExtLimit.lean`の
`exists_descendPieceR_flat_mvPolynomial_baseChange`(この作戦で3回
ビルドして3回とも別の場所で落ち、`FieldLimit.lean`の
`localization_away_quotient_mvPolynomial_flat_equiv`系を`≃ₐ[B']`へ
作り直す方針に切り替えた)。

## 53. `lake build ABC3` は `CorrHyp/` を**ビルドしない**——検証の的を
間違えると「0エラー」が何の証拠にもならない(2026-09-05)

**事実**: `lean/lakefile.toml` の `[[lean_lib]] name = "ABC3"` には
glob 指定が無いので、`lake build ABC3` が作るのは **`ABC3.lean`
(ルート)から `import` で辿れるモジュールだけ**である。ルートは
`Meta.Claim / Meta.Calibration / Interface / Skeleton / Found / Gap /
Check` の7本を import するが、`ABC3/Found.lean` には `CorrHyp` の行が
1つも無く、`ABC3/Found/CorrHyp.lean` という集約モジュールも存在しない。
つまり `Found/CorrHyp/*.lean` は**ルートから到達できない**。

**症状**: `lake build ABC3` が「Build completed successfully (6590
jobs)」と言っても、`CorrHyp/FieldLimit.lean` や `ExtLimit.lean` に
書いた定理は**一度もコンパイルされていない**。それどころか、明示
ターゲットのビルドに失敗した直後は `ExtLimit.olean` が**消えたまま**に
なり、`lake build ABC3` は何事も無かったかのように成功する。この状態で
MCP の `lean_start(["ABC3.Found.CorrHyp.ExtLimit"])` を呼ぶと、
3秒ほどで「読み込んだ」と返るのに中身は空で、
`unknown namespace AlgebraicGeometry` や `Γ(X, U)` の parse error
(`unexpected token '('; expected ')'`)という**一見無関係な**エラーに
なる——「環境が壊れた」ように見えるが、実体は olean が無いだけ。

**正しい検証の的**: `CorrHyp` を触ったら
`lake build ABC3.Found.CorrHyp.Instance4`(この系列で最も深いモジュール、
`ExtLimit → FieldLimit → SchemeFEt / QcqsSpace` を全部含む)を明示的に
叩く。`lake build ABC3` は**それ以外**の部分の回帰検査として併用する。

**関連**: `lake build` 実行中に MCP の `lean_start`/`lean_check` を
並行で走らせると、書き換え途中の olean を読んで同じ「空の環境」症状に
なる。ビルド中は REPL を触らないこと。

実例: 2026-09-05 のセッションで、`ExtLimit.lean` の巨大定理が3回失敗
した直後、`lake build ABC3` は 0 エラーのままだったが `ExtLimit.olean`
は存在せず、REPL が空環境になっていた。

## 54. **「almost split ⟹ almost 消滅」**——`c•𝟙` が消える対象を経由する
と分かれば、コホモロジーの `c` 零化は関手性だけで出る(2026-09-05、
1 セッションで 2 回使えた再利用パターン)

**状況**: 「`m`(あるいは `p^ε`)が `H^*(…)` を零化する」型の主張
(almost mathematics で頻出。Faltings の remark 2.1(v) の Hochschild
cohomology、Theorem 2.4(ii) の群コホモロジーの両方がこの形)。

**やらなくてよいこと**: コチェイン複体の水準で縮約ホモトピーを書く
(`inhomogeneousCochains` の `Fin.contractNth` を相手にする)。低次
(`H^1`・`H^2`)なら明示公式で押せるが、一般次数では相当重い。

**やること**: 対象 `M` が「`c` 倍だけずれた直和因子」であること、
すなわち **`s : M ⟶ N`・`μ : N ⟶ M` で `s ≫ μ = c • 𝟙 M`** を作り、
`N` 側のコホモロジーが消えることを言う。あとは関手性だけ:
`F(s) ≫ F(μ) = F(s ≫ μ) = F(c•𝟙) = c•𝟙` が `F(N) = 0` を経由するので
`c•𝟙 = 0`。**「`c` で消える」は `c•𝟙` が零射になる、と読み替えるのが鍵**。

適用例(どちらも `lean/ABC3/Found/Falt1/`):
- `ext_smul_eq_zero_of_almost_split`(`AlmostEtale.lean`):
  `N := T`(可換環 `T` 自身、`T` 上射影的)、`Ext^{k+1}(T,M)=0` を
  `Ext.eq_zero_of_projective` で。`Ext.mk₀`・`Ext.mk₀_comp_mk₀_assoc`・
  `Ext.smul_comp`・`Ext.mk₀_smul` だけで閉じる。
- `transfer_groupCohomology_smul_eq_zero`(`GaloisTransfer.lean`):
  `N := Coind_1^G(M)`、`H^{n+1}(G,N)=0` を Shapiro の補題
  (`groupCohomology.coindIso`)+ 自明群の消滅
  (`isZero_groupCohomology_succ_of_subsingleton`)で。

**関手が `Linear` インスタンスを持たないとき**(`groupCohomology.functor`
も `HomologicalComplex.homologyFunctor` も持っていない、2026-09-05 実測):
`F(c•𝟙) = c•𝟙` が直接は言えない。そのときは**元の水準に降りて**、
`π`(コサイクル → コホモロジー)の自然性(`groupCohomology.π_map_apply`)と
`cocyclesMap` が `c•𝟙` を `c•` に送ること(`iCocycles` が mono なので
`HomologicalComplex.cyclesMap_i` で確かめられる)を使う。この2つで
「`c • π x = F(c•𝟙)(π x)`」が言えるので、あとは `F(c•𝟙) = 0` を
morphism の水準で示せばよい(`map_id_comp` + `IsZero.eq_zero_of_tgt`)。

**教訓**: almost mathematics の「`m` が零化する」は、ホモロジー代数
としては**射影性/入射性の `c` 倍版**である。対象そのものの分解を
探すより、`s ≫ μ = c • 𝟙` という 1 本の等式に落とすと、あとは
mathlib の既存の関手性 API に載る。

## 55. `Subring` の `def` を挟むと `CommRing` インスタンス経路が変わり `rw` が落ちる(`▸` なら通る)

`adjoinIntegers K x`(`def`、`Subring L`)と `𝒪[L]`(`Valued.integer L`)は
`rfl` で一致するが、`↥` を取ったときの `CommRing` インスタンスの**経路**が違う:

```
↥(adjoinIntegers K x) : CommRing.toCommSemiring.toSemiring
↥𝒪[L]                : (SubsemiringClass.toCommSemiring 𝒪[L]).toSemiring
```

このため両者にまたがる書き換えは `rw` が
「Application type mismatch / not type-correct under `instances` transparency」
で落ちる。**`▸` の項レベルのキャスト**(または `exact`)にすると
defeq で通る。#37(`rw` は syntactic・`exact` は defeq)と同じ根。

実例: `Found/PGC/UnramifiedExtension.lean::isAdic_maximalIdeal_adjoinIntegers`
——基礎体版(`ValuationRingComplete.lean`)は `rw [show ... from rfl, hball n]`
で書けたが、拡大体では `have h1 : ... := ...; exact (hball n) ▸ h1` にする必要が
あった。

あわせて: `Valued.integer.norm_irreducible_pos` などは
`[NontriviallyNormedField K]` を要求する。拡大体 `↥K.carrier⟮x⟯` には
その instance が**無い**ので、`letI := nontriviallyNormedField_adjoin K x`
(`@[implicit_reducible]` 付きの `def`)を先に置く。置き忘れると
インスタンス探索が `whnf` で 200000 heartbeats タイムアウトする
(#47「`maxHeartbeats` の `whnf` タイムアウト」と同じ症状——原因は「構造的に不可能」ではなく「instance が単に無い」)。

## 56. 降下補題を書く前に、**そのファイル自身の在庫**を引く——
`FieldLimit.lean` には既に降下の道具一式がある(2026-09-05)

**やらかし**: `A⊗[ℚ]ℝ = colim_R (A⊗[ℚ]R.1)` のフィルター余極限性から
「ℝレベルでイデアルに属していれば有限段階でも属する」を導く補題を
`exists_fgSubalgebra_mvPolynomial_ideal_mem_descend` として新規に書いた。
ところが同じ主張が同じファイルに
**`exists_mem_ideal_span_range_descend`** としてすでに存在しており、
しかも`exists_mvPolynomial_quotient_ringHom_descend2`の証明中で
実際に使われていた。完全な重複だったので削除した。

**在庫の引き方(`CLAUDE.md`の「在庫」の具体化)**: `CorrHyp` の降下まわりは
`FieldLimit.lean` に集中しているので、まず

```
grep -n "^theorem exists_fg\|^theorem exists_mvPolynomial\|^theorem mem_ideal\|^theorem exists_mem" \
  lean/ABC3/Found/CorrHyp/FieldLimit.lean
```

で**名前の一覧**を眺めるだけでよい(30 行ほど)。`decl-index` を作るより速い。
主な在庫(2026-09-05 時点):

- `exists_fg_subalgebra_tensor_finset` / `_mvPolynomial_finset` /
  `_polynomial_family` …… ℝレベルの有限個の元・多項式を単一の`R`へ降ろす
- `exists_fgSubalgebra_upperBound` / `upperBound2` …… `R`たちの上界
- `exists_mem_ideal_span_range_descend` …… **イデアル所属の降下**
- `mem_ideal_span_range_promote` …… 所属を`R'`へ持ち上げる
- `exists_mvPolynomial_eval_descend` …… `ψ`(生成元の行き先)の降下
- `exists_mvPolynomial_quotient_ringHom_descend` / `descend2` ……
  **環準同型そのものの降下**(片方向。`descend2`はイデアル所属版)
- `exists_mvPolynomial_quotient_ringEquiv_descend` / `descend'` /
  `_specIso_descend` …… 両方向(`ψ`・`ψ'`)からの同型の降下
- `algebraTensorMap_val_injective` / `mvPolynomial_map_algebraTensorMap_val_injective`
  …… `ℚ`が体なので`A⊗R.1 → A⊗ℝ`は単射
- `algebraTensorMap_val_comp_inclusion` / `algebraTensorMap_inclusion_comp_inclusion`
  …… `val ∘ inclusion = val` などの配管

**ついでの発見**: `descend2` は単一の底`A`を仮定しているが、
底が`A → A'`と動く場合へは**無料で一般化できる**——関係式を先に
`Algebra.TensorProduct.map φ (AlgHom.id ℚ R.1)` で`A'`側へ押し出して
から`descend2`を`A := A'`で使えばよい。
`(map (id A') (val R)) ∘ (map φ (id R)) = map φ (val R)`
(`Algebra.TensorProduct.map_comp`+`id_comp`/`comp_id`)がその根拠。

## 57. MCP REPL の `lean_start` に**独立した 2 つのルート**を渡すと mathlib ごと読み込みに失敗する(無言)

`lean_start` に `["ABC3.Found.PGC.UnramifiedExtension",
"ABC3.Found.PGC.PadicLogSurjective"]` のように、互いに import 関係の無い
2 モジュールを渡すと、**起動は成功したと報告される**(「起動して import を
読み込んだ (5.0 秒)」)のに、環境には mathlib すら入っていない。症状は

* `‖x‖` が `expected token`(ノルム記法が無い)
* 与えたモジュールの宣言が `Unknown identifier`

起動時間が通常(11〜12 秒)の半分以下なのが唯一の手がかり。
`lean_reset` しても直らない。

対処: **`lean_start` のモジュールは 1 つにする**(必要なら、両方を import
する薄いファイルを先に作って `lake build` し、それを 1 つだけ渡す)。
mathlib のモジュールを併記するのは問題ない
(`["ABC3.Found.PGC.UnramifiedExtension", "Mathlib.FieldTheory.…"]` は動く)
——壊れるのはプロジェクト側のルートが 2 つ以上のとき。

なお、そもそも olean が未ビルドのモジュールを渡した場合も同じ無言の失敗に
なる(`lake build <module>` を先に通すこと)。

### #57 の★原因判明(同日追記)

上の「2 つのルートを渡すと壊れる」は REPL の制限ではなく、**その 2 つが
同時に import できない**ことの現れだった。実際に `lake build` で同じ 2 つを
import するファイルを作ると:

```
error: import ABC3.Found.PGC.PadicLogMul failed,
  environment already contains 'ABC3.Found.PGC.coeff_pow_eq_zero_of_lt'
  from ABC3.Found.PGC.AdjoinIntegers
```

`Found/PGC/AdjoinIntegers.lean`(`PowerSeries` 版)と
`Found/PGC/PadicLogMul.lean`(`Polynomial` 版)に**同名の定理**があり、
Lubin-Tate 系と p 進対数系の二つの枝が同時に使えなかった。
`PadicLogMul.lean` 側を `coeff_polynomial_pow_eq_zero_of_lt` に改名して解消。

教訓:
* **REPL が無言で壊れたら、まず `lake build` で同じ import を試す**
  ——REPL はエラーメッセージを落とすが `lake build` は出す。
* 同じ namespace に汎用的な名前(`coeff_pow_eq_zero_of_lt` のような)を
  置くときは、枝が合流する日を考えて修飾語を付ける。

## 58. 構造体インスタンスのフィールドは **改行区切り**にする——`,` 区切りが `refine ⟨{ … }, ?_⟩` の中でパースに失敗する

`Derivation` や `AlgHom` を `refine` の中で組み立てるとき、

```lean
refine ⟨{ toFun := f, map_one' := ?_, map_mul' := ?_, … }, ?_⟩
--                    ^^ ここで
-- error: unexpected identifier; expected '}'
```

`,` 区切りだと `toFun := f` の直後で「`}` を期待」と言われる。
**改行区切りにすれば通る**:

```lean
refine ⟨{ toFun := f
          map_one' := ?_
          map_mul' := ?_ }, ?_⟩
```

入れ子の構造体があるときは、内側の `}` を**独立した行**に置き、
外側のフィールドを内側の `{` と同じ列に揃える:

```lean
refine ⟨{
  toLinearMap := {
    toFun := fun b => …
    map_add' := …
    map_smul' := …
  }
  map_one_eq_zero' := ?_
  leibniz' := ?_
}, ?_⟩
```

(2026-09-05、`Found/Falt1/AlmostDerivation.lean` の
`exists_derivation_extension` で 2 回踏んだ。)

## 59. **非可換環での多項式恒等式**は `Polynomial.eval₂RingHom'` で移す

行列環のような非可換 `S` では `Polynomial.aeval` が使えない
(`CommSemiring S` を要求する)。しかし係数と可換な元 `u : S` があれば

```lean
Polynomial.eval₂RingHom' (algebraMap A S) u
  (fun r => Algebra.commute_algebraMap_left r u) : A[X] →+* S
```

が**環準同型**になるので、可換環 `A[X]` で `ring` で証明した恒等式を
`congrArg` でそのまま `S` へ移せる。

```lean
have h := congrArg φ (my_poly_identity s)
simp only [hφ, Polynomial.eval₂RingHom'_apply, map_sub, map_mul, map_pow,
  map_ofNat, Polynomial.eval₂_C, Polynomial.eval₂_X] at h
```

`noncomm_ring` は mathlib の該当ファイルを import していないと使えず、
そもそも仮定を使えないので、この経路のほうが強い。
(2026-09-05、Faltings の "Tripling ε"(almost 冪等元の持ち上げ)
`Found/Falt1/AlmostDeform.lean` で使用。)

## 60. `A = B` を代入した具体例で `Semiring.toModule` と `Algebra.toModule` が衝突する

`IsAlmostEtaleCovering (A := X) (B := X) p` のように**同じ型を 2 回**
渡した具体例を作ると、

```
synthesized type class instance is not definitionally equal to
expression inferred by typing rules,
  synthesized  Semiring.toModule
  inferred     Algebra.toModule
```

が出る(`haveI := (hAE.2.1 : Module.Finite _ _)` のところ)。
`Module ℤ`(`AddCommGroup.toIntModule` 対 `TensorProduct.instModule`)でも
同種の衝突が起きる。

**逃げ方**: 具体例では `A ≠ B` を選ぶ(`B := Fin 2 → A` など)、
底環を `ℤ` 以外にする(`Polynomial ℤ` など)。非空虚性の対照としては
むしろ非自明な例のほうが良いので、これは実質的な損にならない。
(2026-09-05、`Found/Falt1/Section3.lean`・`Section4.lean` の
非空虚性で 2 回踏んだ。)

## 61. `Found/<Track>/` のファイルは **`Found.lean` に import を足さないと既定の `lake build` に入らない**

`lake build`(既定ターゲット `ABC3`)がビルドするのは `ABC3.lean` から
**推移的に到達できるモジュールだけ**である。`lean/ABC3/Found/Falt1/*.lean`
を 9 ファイル書いても、`Found.lean` に import が無ければ既定ビルドは
一度もそれらを触らない——`lake build ABC3.Found.Falt1.X` では通るので
気づきにくい。

**症状**: `lake build` が成功しているのに、MCP REPL を再起動しても
新しく書いた定理が「unknown identifier」になる(olean が古いまま)。

**対策**: 新しいトラックのファイルを作ったら、その場で `Found.lean`
(または対応する集約ファイル)に import を足し、`lake build` の job 数が
増えることを確認する。
(2026-09-05 に判明。Falt1 の 9 ファイルが既定ビルドの外にあり、
それまでのコミットの「lake build 6590 jobs 成功」は Falt1 の検証に
なっていなかった。登録後は 6800 jobs。)

## 58. 構造体の Prop フィールドが `rw` を止める(motive is not type correct)

`PAdicLocalField` のように「型 + その型に依存する Prop の証明」を持つ構造体は、
中身の型を `rw` で書き換えられない:

```lean
-- hx : K.carrier⟮x⟯ = IntermediateField.fixedField H
unfold adjoinField fixedFieldLocalField
rw [hx]
-- ✗ motive is not type correct:
--   adjoinField._proof_4 K x : Module.Finite ℚ_[p] ↥K.carrier⟮x⟯
--   but is expected to have type FiniteDimensional ℚ_[p] ↥_a
```

証明項が抽象化されないため motive が型付かない。**Prop なのに落ちる**のが罠。

**直し方**: その Prop を**明示引数**で取る共通のコンストラクタを 1 つ挟む。

```lean
noncomputable def intermediateLocalField (K : PAdicLocalField p)
    (E : IntermediateField K.carrier K.closure) (hE : FiniteDimensional ℚ_[p] E) :
    PAdicLocalField p := { carrier := E, isFinite := hE }

theorem intermediateLocalField_congr (K) {A B} (hA) (hB) (h : A = B) :
    intermediateLocalField K A hA = intermediateLocalField K B hB := by
  subst h; rfl        -- 証明無関係で hA = hB は defeq
```

元の 2 つの定義がこの形に `rfl` で一致することを 1 行ずつ示せば、
以後は `rw [thisEq, thatEq]; exact intermediateLocalField_congr K _ _ h` で通る。
実例: `Found/PGC/UnramifiedCriterion.lean`(2026-09-05)。

## 62. MCP REPL は `autoImplicit` が効いている——`Type v` が REPL で通ってもファイルで落ちる

プロジェクトの `lakefile.toml` は

```toml
relaxedAutoImplicit = false
autoImplicit = false
```

だが、**MCP REPL の基準環境はこの設定を引き継がない**。そのため

```lean
theorem foo {R : Type*} [Ring R] :
    ∀ (r : ℕ) (N : Fin r → Type v), …   -- ★ v を宣言していない
```

は REPL では `v` が自動束縛されて通るが、ファイルに書いて `lake build`
すると `unknown universe level 'v'` で落ちる。

**対策**: 宇宙変数は `Type*` を使うか、ファイル冒頭で `universe v` を
宣言する。より一般に、`lean_check` で通ったものは**必ず**
`lake build <module>` で確認してからコミットする(#50(b) の一般則の
具体例)。

(2026-09-05、`Found/Falt1/Section1.lean` の `length_pi_fin` で踏んだ。)

## 63. `MvPolynomial.eval₂_comp_left` が「パターンが見つからない」——ゴールが素の `eval₂`、補題は `RingHom` 適用形

`eval₂_comp_left (k : S →+* T) (f) (g) (p) : k (eval₂ f g p) = eval₂ (k.comp f) (k ∘ g) p`
の左辺は **`k` の coe 適用**である。一方、ゴールに現れる外側は
`MvPolynomial.eval₂ algV valV (…)` という**素の関数適用**で書かれていることが多い。
`eval₂ f g = ⇑(eval₂Hom f g)` は `rfl` だが**構文が違う**ので `rw` は噛まない。

```
error: Did not find an occurrence of the pattern
  (MvPolynomial.eval₂Hom algV valV) (MvPolynomial.eval₂ ?f ?g ?p)
in the target expression
  MvPolynomial.eval₂ algV valV (MvPolynomial.eval₂ MvPolynomial.C ψ (...)) = 0
```

**直し方**: `show` で頭を `eval₂Hom` の適用形に揃えてから `rw`。defeq なので通る。

```lean
show (MvPolynomial.eval₂Hom algV valV) (MvPolynomial.eval₂ MvPolynomial.C ψ …) = 0
rw [MvPolynomial.eval₂_comp_left (MvPolynomial.eval₂Hom algV valV)]
```

**同じ回の第2の穴**: 続けて `rw [← eval₂_comp_left restr algU valU p]` が
`(fun i => restr (valU i))` と `(⇑restr ∘ valU)` の食い違いで落ちた。
`funext` で作る補助等式の**型の方を `⇑restr ∘ valU` で書く**と一発で揃う
(項は同じなので `funext hψval` はそのまま通る)。

```lean
have hfun : ((MvPolynomial.eval₂Hom algV valV : MvPolynomial ι' 𝔹 →+* T) ∘ ψ)
    = (⇑restr ∘ valU) := funext hψval   -- 型を ∘ の形で書くのがコツ
```

実例: `Found/CorrHyp/FieldLimit.lean` の `eval₂_map_aeval_eq_zero`。

## 64. テンソル積からの環準同型の等式に `ext x` を使うと `TensorProduct.ext` が先に噛む

`f g : (A ⊗[R] B) →+* T` の `f = g` を示すのに `ext x` と書くと、
`Algebra.TensorProduct.ext`(成分ごとに見る `ext` 補題)が先に選ばれて
`x : A` になり、続く `induction x using TensorProduct.induction_on` が

```
Invalid target: x has type ↑Γ(X.left, U) but is expected to have type ?m ⊗[?m] ?m
```

で落ちる。**直し方**: `ext` を使わず `refine RingHom.ext fun x => ?_` と
明示する。そのあと `induction x using TensorProduct.induction_on` が通る。

**`add` ケースの注意**: 帰納法の仮説 `hx`/`hy` は `(f.comp g) x = (h.comp k) x`
という **`.comp` の形**で出るのに、`simp` はゴールを `f (g x) = h (k x)` へ
展開してしまうため `simp [hx, hy]` が「引数が使われていない」と言って
失敗する。両方を同じ形に揃えてから `rw` する:

```lean
| add x y hx hy =>
    simp only [RingHom.comp_apply, map_add] at hx hy ⊢
    rw [hx, hy]
```

実例: `Found/CorrHyp/ExtLimit.lean` の `pieceAlgebraMap_comp_naturality`
(純テンソル上の等式 `pieceAlgebraMap_naturality` を環準同型の等式へ
持ち上げるところ)。

## 65. 巨大な型では `rw` が `whnf` のヒートビートを食い尽くす——`congrArg` + `Eq.trans` の**項**で組む

`Γ(C, α ⁻¹ᵁ (pullback.fst X.hom toBaseK ⁻¹ᵁ U))` のように `letI` を
何段も重ねた型の上では、`rw` の motive 計算と `whnf` が爆発する。

```
error: (deterministic) timeout at `whnf`, maximum number of heartbeats (4000000) has been reached
```

`maxHeartbeats 4000000` でも足りない。同じ内容を**項で組む**と
`maxHeartbeats 1000000` のまま **2 秒**で通る:

```lean
-- ✗ 87 秒かけて timeout
    rw [← hq k, MvPolynomial.map_map]
    exact congrArg (fun t => MvPolynomial.map t.toRingHom q₀) hsplit

-- ✓ 2 秒
    Eq.trans (congrArg (fun t : A →ₐ[ℚ] B => MvPolynomial.map t.toRingHom q₀) hsplit)
      (Eq.trans (MvPolynomial.map_map _ _ _).symm
        (congrArg (MvPolynomial.map g.toRingHom) (hq k)))
```

`congrArg` の関数には**型注釈を付ける**(`fun t : A →ₐ[ℚ] B => …`)。
付けないと `t` の型が決まらず、結局同じ探索に落ちる。

**併せて**: 補題を述べるとき `algebraMap R S` と書かず
`(Algebra.TensorProduct.map (AlgHom.id ℚ A) (Subalgebra.val R.1)).toRingHom`
のような**明示形**にしておくと、使う側で `algebraMap` の
`RingHom.toAlgebra` 越しの defeq 判定が起きず、これも爆発を防ぐ。
実例: `Found/CorrHyp/ExtLimit.lean` の
`pieceAlgebra_relation_descend_q₀_map` と `descendPieceR_hψ`。

## 59. 中間体の 2 層をまたぐ `rfl` は kernel を止める(1 層なら速い)

`K.closure` の中間体 `K⟮x⟯ ≤ K⟮y⟯` に沿った写像を作るとき:

```lean
-- ✗ IntermediateField.inclusion を使うと kernel deterministic timeout(実測 60 秒)
example (z : ↥K⟮x⟯) : ((IntermediateField.inclusion hle z : ↥K⟮y⟯) : K.closure)
    = ((z : ↥K⟮x⟯) : K.closure) := rfl

-- ✓ 素直に書けば 0.06 秒(中間体の元は 1 層)
example (z : ↥K⟮x⟯) : ‖(⟨(z : K.closure), hle z.2⟩ : ↥K⟮y⟯)‖ = ‖z‖ := rfl
```

ところが `adjoinIntegers K x`(中間体の**部分環**)の元は **2 層**なので、
同じ形の `def` を書いても kernel が落ちる(実測 60 秒)。
`maxHeartbeats 2000000` を付ければ `def` 自体は 30 秒で通るが、
そこから `norm ... = norm ...` を `rfl` で出そうとするとまた落ちる。

**回避策(実測で効いた)**: ノルム保存を**明示補題として先に**用意し、
`def` の中では `rw` でそれを使う——`z.2` を defeq に頼らせない。

```lean
theorem norm_mk_of_le (w : ↥K⟮x⟯) : ‖(⟨(w : K.closure), hle w.2⟩ : ↥K⟮y⟯)‖ = ‖w‖ := rfl  -- 1 層、速い

noncomputable def adjoinIntegersIncl (z : adjoinIntegers K x) : adjoinIntegers K y :=
  ⟨⟨((z : ↥K⟮x⟯) : K.closure), hle _⟩,
   by show ‖(⟨_, _⟩ : ↥K⟮y⟯)‖ ≤ 1
      rw [norm_mk_of_le K hle]      -- ★ここが要
      exact z.2⟩
```

**60 秒 → 0.12 秒**。同じ形で `RingHom` にすると 3.3 秒(defeq に頼る版は
maxHeartbeats 2000000 でも 30 秒)。

実例: `Found/PGC/TotallyRamified.lean` の docstring(2026-09-05)。

## 66. 在庫確認を「結論の形」でやる——`map` と `aeval` の交換律を2度書きかけた

`MvPolynomial.map ψ (aeval ev p) = aeval (map ψ ∘ ev) (map ψ p)` を
新規に証明してから(0.13 秒で通ったので気付かず)、`FieldLimit.lean` に
**同じ主張が `mvPolynomial_map_aeval_comm_general` として既にあった**
ことに気付いた。`mvPolynomial_map_aeval_comm`(`Algebra.TensorProduct.map`
専用版)と、その一般版の 2 本が並んでいた。

**教訓**: `decl-index` を**名前**で引くと当たらない(自分が付けたい名前は
`map_aeval_map`、既存は `mvPolynomial_map_aeval_comm_general`)。
`.cache/decl-index.txt` は statement 付きなので、**結論のリテラル**
(ここなら `aeval` と `map` が両方出てくる行)で grep する:

```
grep -n "^theorem .*aeval" lean/ABC3/Found/CorrHyp/FieldLimit.lean
node tools/decl-index.mjs && grep "aeval.*map\|map.*aeval" .cache/decl-index.txt
```

同じ回に `algebraTensorMap_inclusion_comp_inclusion` も既存だと分かり、
新規に要ったのは「係数側も動く版」`_of_map` 1 本だけだった。

## 67. Git Bash の `python` は Microsoft Store のスタブ——**黙って何もしない**

`python - << 'PYEOF' … PYEOF` でファイルを書き換えたつもりが、**2 回連続で
何も書かれていなかった**。原因は PATH 上の

```
/c/Users/Aruta/AppData/Local/Microsoft/WindowsApps/python
```

が Microsoft Store のインストーラ用スタブで、標準入力を読まずに終了するため。
終了コードも 0 になりうるので `&&` チェーンでも気付けない。

**直し方**:
- ファイル追記は `cat >> file << 'EOF'`(シェル組み込み)を使う。
- Python が要るときは CLAUDE.md に書かれている絶対パス
  `C:\Users\Aruta\miniforge3\envs\py311env\python.exe` を使う
  (`PYTHONIOENCODING=utf-8` も忘れずに)。
- **検算**: `git diff --cached --stat` に意図したファイルが**現れているか**を
  必ず見る。これが今回の発覚の手がかりだった。

同種の「黙って失敗する」罠は #(シェルの終了コード隠蔽)と同じ族——
ログや stat の**実物**を見るまで成功と見なさないこと。

## 68. `Unknown constant` は「mathlib に無い」ではなく「**import していない**」ことが多い

**症状**: `IntermediateField.LinearDisjoint.of_inf_eq_bot` を使ったら
`Unknown constant`。`mcp__abc3-lean__lean_check` の `#check @...` でも同じ。
「pin した mathlib(`db127794`)には無いのだろう」と結論しかけた。

**実際**: `.cache/mathlib-index.txt` を引いたら**あった**
(`FieldTheory/LinearDisjoint.lean:157`)。ABC3 のどこも
`Mathlib.FieldTheory.LinearDisjoint` を import していなかっただけ。
`import Mathlib.FieldTheory.LinearDisjoint` を 1 行足したら通った。

**なぜ嘘に見えるか**: MCP REPL は `Mathlib` 全体ではなく **ABC3 の import 集合**で
動く。だから REPL の `#check` の Unknown も「木に無い」の証拠にならない。
`exact?` が引けないのも同じ理由。

**手順**:

1. `Unknown constant` が出たら、まず
   `grep '<名前>' .cache/mathlib-index.txt`(無ければ
   `node tools/decl-index.mjs --mathlib` で作る)。
2. **あれば** その行のファイルを `import Mathlib.<パス>` で足す(拡張子を除き
   `/` を `.` に)。
3. **無ければ**そのとき初めて「mathlib 不在」と記録する。

★「不在」の測定を誤ると、既にある数学を数百行書き直すことになる。

## 68. 「補題を適用した形」を**独立した statement として書くこと自体**が高い

#65 は「`rw` が `whnf` を食い尽くす」だったが、同じ族のより厄介な現れ方が
ある——**証明ではなく statement の型検査だけで heartbeats が尽きる**。

実測: `descendPieceR_localization_isOpenImmersion` の結論に
`.comp (e : C ≃+* A).toRingHom` を足しただけの定理を書いたところ、
`maxHeartbeats 4000000` で **264 秒かけて timeout**。証明は
`exact isOpenImmersion_specMap_comp_ringEquiv _ _ (既存定理 …)` の 1 行で
済むのに、**その型を書き下ろすと `e` の型(`letI` を何段も抱えている)を
`whnf` する羽目になる**。

**直し方**: 「適用した形」を独立の定理にしない。汎用の移送補題

```lean
theorem isOpenImmersion_specMap_comp_ringEquiv (φ : A →+* B) (e : C ≃+* A)
    (h : IsOpenImmersion (Spec.map (CommRingCat.ofHom φ))) :
    IsOpenImmersion (Spec.map (CommRingCat.ofHom (φ.comp e.toRingHom)))
```

だけを置き、**使う場所(証明の中)で適用する**——そこでは項がすでに文脈に
あるので `whnf` が要らず安い。

**見分け方**: `letI` を 5 段以上抱えた定義(ここでは
`descendPieceRModel_ringEquivQuotientMap`)を**型の中で合成**しようとして
いるなら、この罠に入っている。

## 69. `adjoinField K x` と `adjoinIntegers K x` の境界は「工夫で越える」ものではない(実測)

#59 の続き。`ABC3` には同じ整数環の**二つの表現**がある:

* `adjoinIntegers K x : Subring ↥K⟮x⟯`(手で作った `{y | ‖y‖ ≤ 1}`)
* `𝒪[(adjoinField K x).carrier]`(`Valued` から来る `Valuation.integer`)

`integersEquivAdjoinIntegers` はこの二つを繋ぐ `≃+*` で、**構成は通る**
(`card_residueField_adjoinField` は実際にこれを使っている)。しかし

```lean
theorem tst (w : 𝒪[(adjoinField K x).carrier]) :
    (integersEquivAdjoinIntegers K x w).1 = w.1 := rfl
```

——定義上ほとんど自明(`toFun z := ⟨z.1, _⟩`)なのに——は**通らない**:

| 設定 | 結果 | 時間 |
|---|---|---|
| 既定 | `(kernel) deterministic timeout` | 212 秒 |
| `maxHeartbeats 1000000` | `(deterministic) timeout at whnf` | 126 秒 |

(2026-09-05 実測。`Subtype.ext rfl` でも同じ。)

**原因**: 両辺の型が `↥K⟮x⟯` と `(adjoinField K x).carrier` で異なり、
`𝒪[·]` がスペクトルノルム由来の `Valued` インスタンスなので `whnf` が
展開しきれない。#59 の「1 層なら速い」の**外側**にある。

**対処**: 越えようとしない。**片側に寄せて書き直す**。

### ★訂正(2026-09-05 追記): 寄せる先は `adjoinIntegers` 側

上の表(212 秒 / 126 秒)は実測そのままで正しい。しかし当初ここには
「`𝒪[(adjoinField K x).carrier]` 側で直接書けば通る」と書いていた——
**寄せる側が逆だった**。

`Gal(K(x)/K)` の整数環・剰余体への作用は、`UnramifiedExtension.lean` に
`algEquivIntegers`・`residueAlgEquiv`・`residueGalHom` として
**すでに `adjoinIntegers` 側で完結して存在する**。`exists_frobenius` が
与える性質(`residueAlgEquiv K x σ z = z ^ q`)も `adjoinIntegers` 側。
それを `𝒪[·]` 側で作り直したのが旧 `UnramifiedFrobenius.lean::integersEquivOf`
/ `residueEquivOf` で、**二重化はそこで発生し、橋を渡ろうとした瞬間に
kernel が止まっていた**。橋は渡る必要が無かった。

**判定法**: 越えられない境界に見えたら、「どちらが橋を必要としないか」を
先に数える。**`𝒪[·]` を一度も書かなければこの壁は発生しない。**
旧 `integersEquivOf` / `residueEquivOf` は削除し、
`Found/PGC/UnramifiedFrobenius.lean` は `algEquivIntegers` /
`residueAlgEquiv` だけで書き直した(`𝒪[(adjoinField K x).carrier]` を
一度も書かない)。書き直した結果は
`algEquivIntegers_eq_pow_of_pow_eq_one` まで通っている。

## 70. `linarith` / `norm_num at h` は**文脈全体**を前処理する——巨大な型の項が居ると止まる

`↥K⟮x⟯`(スペクトルノルム由来の `NormedField` インスタンス)の元が
文脈に居るゴールで、最後の 1 行が算数でも次のように書くと**返って来ない**:

```lean
· have h2 : (0:ℝ) ≤ ‖ζ‖ := norm_nonneg ζ
  rw [h.1] at h2      -- h2 : (0:ℝ) ≤ -1
  linarith            -- ← 400 秒経っても返らない(REPL を落とすしかない)
```

`norm_num at h2` でも同じ(600 秒で打ち切り、2026-09-05 実測)。
`h2 : (0:ℝ) ≤ -1` は自明なのに止まるのは、`linarith` / `norm_num at` が
**ゴールだけでなく文脈の全仮説を前処理する**ため——同じ文脈に
`ζ : ↥K⟮x⟯`・`hζ : ζ ^ m = 1`・`h1 : ‖ζ‖ ^ m = 1` が居ると、
そこで `whnf` が爆発する。

**対処**: 文脈を見ないで済む形、つまり**項**で書く。同じ枝が **0.82 秒**になる:

```lean
· exact absurd (h.1 ▸ norm_nonneg ζ) (by norm_num)
```

`by norm_num` のゴールは `¬ ((0:ℝ) ≤ -1)` だけで、`at` が付いていないので
文脈を触らない。#65(巨大な型では `rw` が heartbeats を食う)の
「自動化タクティクは文脈のサイズに比例して重い」という同じ話の、
決定手続き側の現れ。

**目印**: `Found/PGC/` の `↥K⟮x⟯`・`K.closure`・`𝒪[·]` が文脈に居る証明で
`linarith` / `nlinarith` / `positivity` / `norm_num at` / `omega` を
書こうとしたら、まず項で書けないかを見る。

## 71. 同型な2つの構造を「別の項」として作るには、**台の型は同じまま**にして型同義語で包む

「同じ型の上に別の代数構造を2つ載せ、その2項が異なることを証明したい」ときの型。
`Check/PGC/Prop12Degenerate.lean`(Prop 1.2 の反証)で踏んだ。

**やってはいけない**: `ULift ℚ_[p]` のような **`structure` で包んだ別の型**にすること。
`ULift.field` は在る(2026-08-14 の「mathlib に無い」は今は誤り)が、
`ULift ℚ_[p] ≠ ℚ_[p]`(**型の非等号**)は Lean では証明できないので、
「2つの項が異なる」が言えない。

**やってはいけない(その2)**: 台を `ℚ_[p]` のまま、捻った `Field`/`Algebra` を
その場で使うこと。台の型が標準構造の型と同じなので、`AlgEquiv.symm` や
`Module.Finite.equiv` のインスタンス引数を**探索が標準構造で埋めてしまい**、

```
error: synthesized type class instance is not definitionally equal to
expression inferred by typing rules
```

で落ちる(単一化が割り当てた捻り側と、探索が見つけた標準側が食い違う)。

**正解**: `def TwistedQp (p) : Type := ℚ_[p]`(**`abbrev` ではなく `def`**)。

* `TwistedQp p = ℚ_[p]` は `rfl` で証明できる(項の区別に使える)
* インスタンス探索の鍵は `TwistedQp` なので、捻ったインスタンスだけが見つかる
  ——`.symm` も `Module.Finite.equiv` も素直に通る

★ただし**捻ったインスタンスを作る所だけ**は台が `ℚ_[p]` のままなので、
`Equiv.algebra` / `Equiv.algEquiv` の引数を**全部 `@` で明示**する必要がある
(`by infer_instance` を書く)。明示せずに型同義語側で直接作ると、
`whnf` が 200000 heartbeats を食って `Type mismatch` になる。
順序は「① 台 `ℚ_[p]` の上で `@` 明示で作る → ② 型同義語のインスタンスとして
`:= negField p` のように**一段の delta で**渡す」。

## `n` を `n-1+1` に `rw` すると、台の型が `n` に依存していて motive が壊れる（2026-09-05）

`ℤ_[p]` の中で二項展開の和 `∑ k ∈ range (p+1), f k` から先頭2項を剥がしたい。
`Finset.sum_range_succ'` を2回使うため `p = (p-1)+1` を `rw` したくなるが、

```
rw [hpe]   -- hpe : p = p - 1 + 1
-- ✗ motive is not type correct:
--   Application type mismatch: hp : Fact (Nat.Prime p) but expected Fact (Nat.Prime _a)
--   in the application @PadicInt _a
```

`p` は `ℤ_[p]` と `Fact p.Prime` の**両方**に現れるので、`p` を抽象化した motive に型が付かない。

**正解**: 剥がす操作だけを**台に依存しない補題**として切り出し、そこで `obtain ⟨m, rfl⟩` する。

```lean
private theorem sum_range_split_two {M : Type*} [AddCommMonoid M] {n : ℕ} (hn : 2 ≤ n)
    (f : ℕ → M) :
    ∑ k ∈ range (n+1), f k = (∑ k ∈ range (n-1), f (k+2)) + f 1 + f 0 := by
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 2 := ⟨n - 2, by omega⟩   -- ★ここでは n に依存する型が無い
  rw [Finset.sum_range_succ' f (m+2)]; congr 1
  rw [Finset.sum_range_succ' (fun k => f (k+1)) (m+1)]; simp
```

呼ぶ側は `rw [sum_range_split_two (by omega) f] at h` だけで済み、`p` は一切書き換わらない。
★一般に「添字の算術を直す `rw`」は、台が添字に依存していると必ずこれを踏む。
**添字の算術は台を `Type*` に一般化した補題の中に閉じ込める。**

## `(selfField p).carrier` は `exact` では `ℚ_[p]` と合うが `rw` では合わない（2026-09-05）

`selfField p : PAdicLocalField p := { carrier := ℚ_[p] }` は `def`（semireducible）なので、
`(selfField p).carrier` は **default 透明度でだけ** `ℚ_[p]` に落ちる。

* `exact`／項の適用は通る —— `ℚ_[p]` について述べた定理を
  `y : (selfField p).carrier` にそのまま当てられる（`tst_qp_odd p hp3 h` で通った）
* `rw` は落ちる —— インスタンス引数を `instances` 透明度で照合するため
  「`Membership (selfField p).closure (Set (AlgebraicClosure ℚ_[p]))`」のような
  型不一致になる（§1 と同じ症状）

**書き方**: 補題は `ℚ_[p]` 側で述べ、`selfField` 側の証明では `rw` を使わず
`exact` / `apply` で当てる。`rw` が要るなら、両辺とも `(selfField p).carrier` 側の
語彙（`IntermediateField.mem_bot (F := (selfField p).carrier)` のように `F`/`E` を明示）で書く。

## 全体ビルドの `failed to read file '….olean.private'` は REPL のメモリ圧（2026-09-05）

`lake build ABC3`（6900 ジョブ）の途中で

```
error: ABC3/Interface.lean:1:0: failed to read file
  '….lake/build/lib/lean/Mathlib/RingTheory/Jacobson/Ring.olean.private'
```

が出た。**Lean のエラーではなく Windows の mmap 失敗**である——
MCP の `repl.exe` が 3 GB 保持したまま常駐していて、olean の mmap が取れなくなっていた。

**直し方**: 全体ビルドの前に `taskkill //F //IM repl.exe`（`lean_check` は
`lean_start` で建て直せる）。再実行したら同じ木がそのまま通った（6921 ジョブ成功）。
★olean を消したり `lake clean` したりする前に、まずメモリを疑うこと。

## `variable (p : ℕ)` が明示のファイルで `def foo (X : T p)` を作ると `foo p X`（2026-09-05）

`Check/PGC/Prop12Degenerate.lean` は `variable (p : ℕ) [Fact p.Prime]`（**明示**）。
そこに `def residueCardAndDegreeObjectOld (RD : ResidueCardinalityOld p) : …` を足すと、
`p` は自動束縛されて**第1引数**になる。呼ぶ側で `residueCardAndDegreeObjectOld RD` と
書くと

```
Application type mismatch: RD has type ResidueCardinalityOld p … but is expected to have type ℕ
```

`Skeleton/` 側の対応物（`variable {p : ℕ}` = 暗黙）からコピーしてくると必ず踏む。
**ファイルの `variable` の明示／暗黙を先に見る。**

## `set A := adjoinIntegers K x` は `IsLocalRing`/`ResidueField` のインスタンス探索を壊す（2026-09-05）

**症状**: `set A := adjoinIntegers K x with hA` としてから
`IsLocalRing.ResidueField A` を使うと

```
Application type mismatch: the argument instIsLocalRingAdjoinIntegers K x has type
  @IsLocalRing (↥(adjoinIntegers K x)) (SubsemiringClass.toCommSemiring ...).toSemiring
but is expected to have type
  @IsLocalRing (↥A) CommRing.toCommSemiring.toSemiring
```

`set` が `↥(adjoinIntegers K x)` を `↥A` に置き換えた結果、`Subring` 由来の
`CommRing` インスタンスが `A.toCommRing`（局所定義の方）として拾われ、
既存の `IsLocalRing` インスタンスと**同じ型に見えなくなる**。副作用として
`residue A a = residue A b` 型の仮説が `a✝` を含む別物になり、`obtain` の
結果が使えなくなる。

**直し方**: 部分環・部分体を `set` で略記しない。`adjoinIntegers K x` を
そのまま書く（長いが 0.2 秒で通る）。略記したいのは剰余体の**位数**の方
なので、`set Q := Fintype.card (IsLocalRing.ResidueField (adjoinIntegers K x))`
のように**ℕ の値だけ**を `set` する。

## `Nat.dvd_sub'` は無くなった（2026-09-05）

`Nat.dvd_sub' : k ∣ m → k ∣ n → k ∣ m - n` は現行 mathlib に無い。
`Nat.dvd_sub`（同じ結論、ℕ の切り捨て減算のまま）を使う。
`p ∣ Q → ¬ p ∣ (Q-1)` は `Nat.dvd_sub hpQ hcon` で `p ∣ Q - (Q-1)` を作り、
`rw [show Q - (Q - 1) = 1 by omega]` で `p ∣ 1` に落とすのが最短。

## `omit [Inst] in` は docstring より**前**に置く（2026-09-05）

```lean
/-- ... -/
omit [CharZero F] in
lemma foo ...        -- ✗ `unexpected token 'omit'; expected 'lemma'`

omit [CharZero F] in
/-- ... -/
lemma foo ...        -- ✓
```

section variable の instance 引数（`[Field F]`・`[CharZero F]` など）は
「`F` を使った宣言」に**自動で全部入る**。1 つでも使わない宣言があると
`unusedSectionVars` 警告が出るので、`omit ... in` を足すか、
`[CharZero F]` を必要な宣言だけの `section` に分ける。

## 型の**添字**に現れる `n` を `rw` しない（2026-09-05）

`f : Γ →* ↥(rootsOfUnity n Ω)` が文脈にあるとき、ゴールの `n` を
`rw [← card_rootsOfUnity hζ hn]`（`n = Nat.card ↥(rootsOfUnity n Ω)`）で
書き換えると `motive is not type correct`（`f` の型の `n` まで抽象化される）。

**直し方**: ゴールではなく**仮説の側**を書き換える。
```lean
have h2 : Nat.card ↥f.range ∣ Nat.card ↥(rootsOfUnity n Ω) := ...
rw [card_rootsOfUnity hζ hn] at h2   -- ✓ h2 : ... ∣ n
```
同じ理由で `IsPrimitiveRoot.pow` に渡す `n = k * d` も
`rw [hk, mul_comm]` ではなく `hk.trans (mul_comm _ _)` と**項で**書く。

## 名前が違う 3 つ（2026-09-05、pGC 経路 C で実測）

| 探した名前 | 実在する名前 |
|---|---|
| `Subgroup.equivOfEq (h : H = K) : H ≃* K` | `MulEquiv.subgroupCongr` |
| `Subgroup.eq_of_le_of_card_le` | `Subgroup.eq_of_le_of_card_ge`（`H ≤ K → card K ≤ card H → H = K`） |
| `AddEquiv.toMultiplicative' : (α ≃+ Additive β) ≃ (Multiplicative α ≃* β)` | `AddEquiv.toMultiplicativeLeft` |

## `ZMod n` に `TopologicalSpace` は無い（2026-09-05）

`ContinuousMonoidHom G (Multiplicative (ZMod n))` は**書けない**。
「連続」を `IsOpen (MonoidHom.ker f)` で定義すれば係数群の位相を一切参照せず
に済む（`Found/PGC/ContinuousHomCount.lean` の `contHom`）。
`ker (f*g) ⊇ ker f ⊓ ker g` から開性を出すには `Subgroup.isOpen_mono` が要り、
これは `[SeparatelyContinuousMul G]` を要求する（`IsTopologicalGroup` から出る）。

## 暗黙の `{m : ℕ}` が `by` ブロック側から**逆に**決まる（2026-09-05、pGC G1 で実測）

`exists_orderOf_eq_of_dvd_card {m : ℕ} (hdvd : m ∣ Nat.card G)` に
`(by rw [card_primeToPTorsion]; exact hdvd)` を渡すと、`m` は `n` ではなく
`Nat.card G` に決まってしまう（`rw` が `m ∣ Nat.card G` を `rfl` で閉じ、
`exact hdvd` が `No goals to be solved` になる）。証明項が `by` ブロックだと
`m` がメタ変数のまま残り、`Nat.card G ∣ Nat.card G` で先に埋まるため。

**直し方**: 暗黙引数を**名前で先に固定**する。
```lean
exists_orderOf_eq_of_dvd_card (G := primeToPTorsion K) (m := n)
  (by rw [card_primeToPTorsion]; exact hdvd)   -- ✓
```

## `Skeleton/PGC/Setup.lean` の `closure` / `absGal` は `abbrev`（2026-09-05）

`K.closure = AlgebraicClosure K.carrier`、`K.absGal = K.closure ≃ₐ[K.carrier] K.closure`
はどちらも `abbrev`（reducible）。`AlgebraicClosure F ≃ₐ[F] AlgebraicClosure F` に
ついて述べた補題は `K.absGal` に**そのまま当たり**、
`contHomCard K.absGal n = contHomCard (AlgebraicClosure K.carrier ≃ₐ[K.carrier] _) n` は
`rfl` で通る（位相・群の instance も一致）。★`autCongrContinuousMulEquiv` /
`contHomCard_congr` で移送する設計を**先に考えない**こと。

## 有限次中間体から**生成元の有限集合**を取る（2026-09-05、pGC D3 で実測）

`FiniteDimensional K ↥E'` から「`adjoin K S = E'` なる有限集合 `S`」を取る道は
`IntermediateField.fg_of_finiteDimensional` **ではない**（そんな名前は無い）。
`fg_of_noetherian` は `[IsNoetherian F E]`（**大きい方の体**全体）を要求するので
中間体には当たらない。実際に効くのは `EssFiniteType` 経由：

```lean
obtain ⟨S, hSfin, hSadj⟩ := IntermediateField.fg_def.mp
  (IntermediateField.essFiniteType_iff.mp inferInstance : E'.FG)
haveI : Finite S := hSfin
```
`Module.Finite → Algebra.FiniteType → Algebra.EssFiniteType.of_finiteType` が
instance で繋がるので `inferInstance` で通る。

★これで「Krull 位相の近傍を**生成元**に落とす」書き方ができる。
`(adjoin F S).fixingSubgroup = fixingSubgroup (Ω ≃ₐ[F] Ω) S`
（`Found/PGC/AdjoinFieldClosure.lean::fixingSubgroup_adjoin_eq`、
証明は `← IntermediateField.le_iff_le` と `adjoin_le_iff` の 2 手）を挟むと、
底体 `F` を `E` に取り替えても条件が「`S` を各点固定する」のまま変わらないので、
**合成体 `E ⊔ E''` も塔も要らなくなる**（`finiteDimensional_sup` は両方の
有限性を要求するので無限次では使えない。そこを回避できる）。

★`_root_.mem_fixingSubgroup_iff` は `M` が**明示引数**：
`(_root_.mem_fixingSubgroup_iff _).mp hf` と書く。`_root_.` を落とすと
`IntermediateField.mem_fixingSubgroup_iff` に取られ、`_` を落とすと
`Unknown constant mem_fixingSubgroup_iff.mp` という**紛らわしいエラー**になる。

## MCP `lean_check` は通るのに `lake build` が `unknown tactic`(2026-09-05)

MCP の基準環境は mathlib を丸ごと読んでいるので `group` / `ring` / `nlinarith` などが
いつでも使える。ところが**自分のファイルの import が狭い**と、ディスクに書いて
`lake build` した瞬間に

```
error: unknown tactic
error: unsolved goals      ← 直後にこれが続く(タクティクが無いので何もしていない)
```

になる。`Found/` 直下に「mathlib だけを import する一般補題」を置くときに必ず踏む。
→ 使ったタクティクの分だけ `import Mathlib.Tactic.Group` / `Mathlib.Tactic.Ring` などを
足す。★`Unknown constant`(#68)の**タクティク版**だが、症状が
「名前が無い」ではなく「unknown tactic + 直後の unsolved goals」なので気づきにくい。
逆に言えば、`lean_check` が通った証明は import を足すだけで必ず通る。
(2026-09-05、`Found/HerbrandIndex.lean`(Herbrand 商)で実測。)

## `(fixedFieldLocalField K H hH).carrier` にはインスタンスが付かない（2026-09-05）

`(selfField p).carrier`（上）の**インスタンス合成版**。`fixedFieldLocalField` は `def`
（semireducible）なので、TC 解決（`instances` 透明度）は射影を潰せない。

```
failed to synthesize instance of type class
  Module K.carrier (fixedFieldLocalField K H hH).carrier
```

`Algebra ℚ_[p] F.carrier` は `PAdicLocalField.isAlgebra` が `?K.carrier` に直接
マッチするので通るが、**底体が `K.carrier` の側**（`Module K.carrier ↥(fixedField H)`）は
射影を剥がさないと見えない。

**直し方**: 補題は生の型 `↥(IntermediateField.fixedField H)` で述べ、`show` で橋を架ける
（`Module.finrank ℚ_[p] (fixedFieldLocalField K H hH).carrier
= Module.finrank ℚ_[p] ↥(IntermediateField.fixedField H)` は `rfl` で通る）。
`Found/PGC/DegreeTransport.lean::finrank_fixedFieldLocalField` が実例。

## `WithZero (Multiplicative ℤ)` の付値は `WithZero.exp` / `WithZero.log` で書く（2026-09-05）

mathlib の adic valuation の戻り値が `Multiplicative.ofAdd (-1)` から
**`WithZero.exp (-1 : ℤ)`** に変わっている（`exp a = coe (Multiplicative.ofAdd a)` なので
定義は同じだが、**simp 補題が `exp`/`log` 側にしか付いていない**）。

- `IsDedekindDomain.HeightOneSpectrum.valuation_exists_uniformizer K v : ∃ π, v.valuation K π = WithZero.exp (-1)`
- `valuation_exists_uniformizer'` は **`π : R`**（整の側）を返す。茂・正則関数の側に錨を打つときはこちら。

★`ord := -(WithZero.unzero h).toAdd` のような定義を扱うなら、まず
`ord f = -WithZero.log (v f)`（**`if` 無し**、`log 0 = 0` なので `f = 0` も含めて成立）
を 1 本立てると、以降が `log_exp` / `log_zpow` / `log_mul` で全部落ちる。
場合分けを毎回書くのと差が大きい。
`Found/Divisor/SchemeWeilOrd.lean::ordPt_eq_neg_log` が実例。

## `lean_start` がツール一覧に無い agent での逃げ道（2026-09-05）

sub-agent によっては MCP が `lean_check` だけを見せていて、`lean_check` は
「まだ lean_start を呼んでいない」と返して起動できない。このときは
**スクラッチパッドに小さい `.lean` を書いて `lake env lean <そのファイル>`**（ガード R1 を
コマンドに `#full-check` を含めて抜ける）で代用する。import が 1 本なら **8〜10 秒**で戻るので、
対象ファイル本体を毎回 `lake build` するよりずっと安い。
★対象ファイル自体を import せず、**head -N で先頭を切り出して新ブロックを挿した写し**を
作ると、もとの `variable` の効き方まで含めて同じ 8 秒で検査できる。

## 構造体フィールドが `[inst]` を含むとき、`f _ _ := 0` の**下線の数**を数え違えると別の場所でエラーになる（2026-09-05）

`ord : ∀ (X : Scheme) [IsIntegral X], PrimeDivisorPt X → X.functionField → ℤ` のような
フィールドを `where` で埋めるとき、**instance-implicit 束縛子も下線 1 個を消費する**
（自動挿入されない）。ところが束縛子を減らしすぎても `Pi.instOfNat` が効いて
「関数への `0`」として**通ってしまう**ことがあり、エラーは
`OfNat (IsNormalScheme x → … → WeilDiv x) 0` のように**別のフィールドで**出る。
→ フィールドの `∀` を数え、`[…]` も 1 個ずつ数えて下線を合わせる。
`Check/FrdI/Ex61OrdDegenerate.lean::zeroWeilOrd`（`ord _ _ _ _ := 0`、`div _ _ _ _ := 0`）が実例。

## Python で Lean ファイルを書き換えるときは、**encode してから開く**（2026-09-05）

`𝓞`（U+1D4DE、BMP 外）を Python のリテラルで `\ud835\udcde` と書くと
**サロゲート対**になり、`utf-8` で encode できない。ところが

    io.open(path, "w", encoding="utf-8").write(out)   # ← これが危険

は **開いた時点でファイルを 0 バイトにし**、そのあと `write` で
`UnicodeEncodeError: surrogates not allowed` を投げる。結果は
**`Found/` のファイルが中身を失ったまま残る**
（実測 2026-09-05: `Found/GenEll/VeluSemistableJ.lean` 117 行 → 0 行。`git checkout --` で復旧）。

→ 直し方は 2 つ。リテラルは **`\U0001D4DE`**（大文字 U の 8 桁）で書くか、
変数に入れて `.replace("𝓞", O)` で差し込む。そのうえで

    data = out.encode("utf-8")        # ★encode してから初めて開く
    with open(path, "wb") as f: f.write(data)

の順にすれば、失敗してもファイルは無傷である。
★先頭に `assert len(src) > 1000` を置いておくと、一度空にしたファイルを
再度上書きして履歴ごと失うのを止められる。
★★heredoc の中の `\\u` は 1 本に潰れて届くことがあるので、
文中にエスケープを書きたいときは `BS = chr(92)` を組み立てて差し込む。

## agent からの 1 循環は `tools/leanfile.mjs` が一番安い（2026-09-05）

`lean_start` がツール一覧に無い sub-agent では、スクラッチパッドに `.lean` を書いて

    node tools/leanfile.mjs <スクラッチの絶対パス>

とする。`lake env lean` を直叩くのと違い `lakefile.toml` の `[leanOptions]`
（`autoImplicit=false` 等）を渡すので `lake build` と食い違わず、**olean を書かない**ので
並行セッションと同じワークツリーでも安全である。
実測: `Found/GaloisRep/Lemma35Ineq` と `Found/GenEll/VeluSemistableJ` を import した
スクラッチ 1 枚で **11〜13 秒**。同じ内容を `lake build` で確かめると
**6 分 45 秒**（4539 ジョブ）だった。
★ただし olean を読むので、**先に直したファイルの新しい宣言**は `lake build` するまで見えない。
新宣言を使う側を検査するときは、両方を 1 枚のスクラッチに並べて書く。

## `haveI : Field ↥E := inferInstance` を先に置くと、直後の `Algebra ↥E Ω` が見つからなくなる（2026-09-06、pGC ノード F で実測）

```lean
example (K : PAdicLocalField p) : True := by
  haveI : Field ↥(unramifiedClosure K) := inferInstance      -- 無害に見える
  haveI : Algebra ↥(unramifiedClosure K) K.closure := inferInstance
  -- failed to synthesize  Algebra (↥(unramifiedClosure K)) K.closure
```

`haveI` はローカル文脈に**新しい fvar のインスタンス**を積む。`IntermediateField.toAlgebra`
は `E` の体構造(標準のもの)に合わせて `Algebra ↥E L` を作るので、
探索がローカルの `Field` を先に拾うと合わなくなる。2 つ目以降でも同じことが起きる
（`haveI : Algebra …` の直後に `Algebra.IsAlgebraic …` が落ちる）。

★対処: **確認のための `haveI` を並べない**。1 宣言に 1 つずつ `example` で試す。
本番でも「すでにインスタンスがあるもの」を `haveI` で置き直さない。
`open scoped NormedField Valued` のせいだと誤診しやすい（実際は無関係）。

## `¬ Countable NNReal` を示すとき、局所仮説の `Countable NNReal` が拾われない（2026-09-06、Check/FrdI/Ex63DegDegenerate で実測）

```lean
theorem not_countable_nnreal : ¬ Countable NNReal := by
  intro h                      -- h : Countable NNReal（局所インスタンスのはず）
  have hinj : Function.Injective (fun x : ℝ => Real.toNNReal (Real.exp x)) := …
  have : Countable ℝ := hinj.countable
  -- failed to synthesize  Countable { r // 0 ≤ r }
```

`NNReal` は `def NNReal := {r : ℝ // 0 ≤ r}` なので、`Function.Injective.countable` の
`[Countable β]` を解こうとしたエラボレータが `β` を**subtype に unfold した形**で
探索し、局所の `h : Countable NNReal` と syntactic に合わなくなる。
`⟨Real.exp x, _⟩` という anonymous constructor で関数を書くと、`β` が最初から
`{r // 0 ≤ r}` に推論されるので同じ穴に落ちる。

★対処: **インスタンス引数を明示で渡す**。

```lean
  haveI : Countable ℝ := @Function.Injective.countable ℝ NNReal h _ hinj
  exact Cardinal.not_countable_real Set.countable_univ
```

関数側も `Real.toNNReal (Real.exp x)` と書いて `β = NNReal` を固定する。
`Uncountable ℝ` のインスタンスと `Cardinal.not_countable_real` は
`Mathlib.Analysis.Real.Cardinality` にある（`Mathlib.Data.Real.Cardinality` は**無い**）。

## MCP `lean_start` の基準環境は sub-agent 間で共有される——他の agent が import を差し替える（2026-09-06 実測）

`lean_start(["Mathlib.NumberTheory.NumberField.ProductFormula", …])` で立てた環境で
数回 `lean_check` したあと、突然 `open NumberField` が `unknown namespace` を返し始めた。
`lean_status` を見ると `imports:` が **こちらが指定していない** `ABC3.Found.PGC.InertiaKummer`
に変わっていた——並行する別 agent が同じ MCP サーバに `lean_start` を打っている。

★対処: **並行セッションだと分かっているときは最初から `node tools/leanfile.mjs` を使う**。
1 往復 11〜13 秒だが、環境を横取りされない。`lean_check` が急に基本語彙を見失ったら
バグを疑う前に `lean_status` の `imports:` を見る。

## `rw [← ker_π]` が motive not type correct になる——書き換え先が**別の宣言の型**に現れている（2026-09-06、Found/PGC/InertiaKummerBound で実測）

```lean
-- ker_restrictUnramified : (restrictUnramified K).ker = (unramifiedClosure K).fixingSubgroup
have hkerle : ∀ f : ↥((restrictInertia K n).ker), (restrictUnramified K).ker ≤ (↑↑f).ker := …
-- 目標: (unramifiedClosure K).fixingSubgroup ≤ (↑↑f).ker
rw [← ker_restrictUnramified]     -- motive is not type correct
```

`f` の型 `↥((restrictInertia K n).ker)` の中に `(unramifiedClosure K).fixingSubgroup` が
**型として**入っている（`restrictInertia` の余域が `contHom ↥(…fixingSubgroup) _`）ので、
その出現まで抽象化しようとして motive が壊れる。`occs` で絞るより

★対処: **書き換えずに、包含を要素ごとに適用する**。

```lean
refine isOpen_ker_of_factors K f.1.2 (fun σ hσ => hkerle f ?_) (hgf f)
rw [ker_restrictUnramified]       -- 目標が `σ ∈ (restrictUnramified K).ker` なら依存が無く通る
exact hσ
```

同じファイルで踏んだ 2 つ目:`congrArg Subtype.val hab`（`hab` が `(fun f => …) a = (fun f => …) b`）は
**β 簡約されない**まま `h1` に入るので、後の `rw [h1]` が「パターンが見つからない」で落ちる。
`have h1 : <明示的な型> := congrArg Subtype.val hab` と型を書けば β 簡約された形で入る。

## `NumberField ↥(⊥ : IntermediateField ℚ ℂ)` は `infer_instance` で出ないが `⟨⟩` で出る（2026-09-06、Check/GenEll/SSCurveNonvacuous で実測）

```lean
example : NumberField (⊥ : IntermediateField ℚ ℂ) := by infer_instance   -- failed to synthesize
example : NumberField (⊥ : IntermediateField ℚ ℂ) := ⟨⟩                  -- 通る
```

`NumberField` は `[to_charZero]` `[to_finiteDimensional]` を**インスタンス暗黙のフィールド**に
持つ `Prop` クラスなので、クラス探索の対象になっていない。中身の 2 つ
（`CharZero ↥⊥`・`FiniteDimensional ℚ ↥⊥`）は `infer_instance` で出るので、
**匿名コンストラクタ `⟨⟩` を書けばよい**。★数体の witness を作るときに毎回踏む。

同じ場面の 2 つ目: `abbrev Kb : IntermediateField ℚ ℂ := ⊥` は
`Complex.instField` が noncomputable なので **`noncomputable abbrev`** と書く要がある
（`def` だけ noncomputable にしても、`abbrev` の側で落ちる）。

## 部分体の中の数値リテラルを `ℂ` へ落とすのは `norm_num`／`push_cast`／`simp` が**全部無力**（2026-09-06 実測）

```lean
example : ((4096 : ↥(⊥ : IntermediateField ℚ ℂ)) : ℂ) = 4096 := by norm_num   -- ⊢ ↑4096 = 4096 が残る
example : … := by push_cast; ring                                            -- 同上
example : … := by simp                                                       -- `simp` made no progress
```

`SubringClass` には `coe_ofNat` に当たる simp 補題が無い（`natCast_mem`・`ofNat_mem` は在る）。
★対処: **包含を環準同型として書いて `map_ofNat` を当てる**。

```lean
example : ((4096 : ↥K) : ℂ) = 4096 := map_ofNat (K.val.toRingHom) 4096
example : (((11 : 𝓞 L)) : L) = 11 := map_ofNat (algebraMap (𝓞 L) L) 11
```

`↑x` と `K.val.toRingHom x` は defeq なので `exact` で通る（`rw` は形が違うので通らない）。
`push_cast` は `neg`・`div` までは押せるので、残った数値リテラルだけこの手で潰すのが速い。

## インスタンス探索は `def` の射影を**開かない**——`E.W.baseChange M` の `Algebra` は `E.fld` で書く（2026-09-06 実測）

```lean
theorem foo (M : Type) [Field M] [Algebra fldQ M] :
    ∃ Q : (ssCurve11a3.W.baseChange M).toAffine.Point, addOrderOf Q = 5 := …
-- failed to synthesize  Algebra (↥ssCurve11a3.K) M
```

`ssCurve11a3 : SSCurve` は普通の `def` なので、`ssCurve11a3.K` を `fldQ` へ潰すには
delta 簡約が要る。**単一化（`exact` の型合わせ）は開くが、インスタンス探索は開かない。**
★対処: 束縛子の側を射影で書く（`[Algebra ssCurve11a3.fld M]`）。
補題を渡す側も `model11a3_baseChange ssCurve11a3.fld M` と射影で書けば、
残りは単一化なので通る。

## 曲線の等式をまたいで `Point` を運ぶには、曲線を**変数に置いてから `subst`**（2026-09-06）

`W₁ = W₂` は命題の等式なので `W₁.toAffine.Point` と `W₂.toAffine.Point` の間の
輸送が要る。`▸` を直接当てると `addOrderOf (h ▸ Q) = 5` が動かせない。
★対処: 曲線を**全称量化した変数**にして `subst` する。

```lean
theorem exists_point_addOrderOf_five (W : WeierstrassCurve L) (hW : W = model11a3 L) :
    ∃ Q : W.toAffine.Point, addOrderOf Q = 5 := by
  subst hW; exact ⟨model11a3P L, model11a3_addOrderOf L⟩
```

`hW : ssCurve11a3.W.baseChange M = model11a3 M` を渡すだけで
`(ssCurve11a3.W.baseChange M).toAffine.Point` の元が出る。

## `letI` を `∃` の**本体**に書いた命題は、`refine` のあとインスタンス探索が届かない（2026-09-06）

`VeluQuotOK` のように `letI : DecidableEq (M : Type) := …` を項の中に書いた命題を
証明するとき、`refine ⟨⊥, inferInstance, ?_⟩` のあとの目標に残った `let` は
**局所文脈に入らない**ので `failed to synthesize DecidableEq ↥⊥` になる。
★対処: `@` で明示的に渡す。

```lean
exact @exists_point_addOrderOf_five _ _ (fun a b => Classical.propDecidable (a = b)) _ h
```

## `AssociatedObject.transport` のような**構造体フィールドの適用**は `simp`/`simpa` が開かない（2026-09-06、Found/PGC/CyclotomicRecovery で実測）

`Skeleton/PGC/Setup.lean` の `AssociatedObject` は `Obj`/`obj`/`transport` を
フィールドに持つ構造体で、`cyclotomicCharacterObject` は無名コンストラクタで定義されている。
`(cyclotomicCharacterObject).transport α f g'` を `f (α.toMulEquiv.symm g')` に
簡約したいとき、`simpa using this` は**射影を開かず**

```
but is expected to have type
  cyclotomicCharacterObject.transport α (cyclotomicCharacterObject.obj K) g' = …
```

で止まる（`simp [cyclotomicCharacterObject]` も、`where` 記法で作った項では効かないことがある）。

★対処: `show` でベータ簡約後の形を書く。`show` は defeq で通るので確実。

```lean
funext g'
show cyclotomicCharacter K.closure p (α.toMulEquiv.symm g').toRingEquiv
  = cyclotomicCharacter K'.closure p g'.toRingEquiv
```

同じ理由で、`congrFun h (α g)` の結果は `have h2 : … := congrFun h (α g)` と
**型を明示して**受けてから `show` で開くのが早い。


## `(M : Type) × (M : Type)` は**直積ではなく Σ 型**にパースされる（2026-09-06、Found/GenEll/VeluDualJ で実測）

`M : IntermediateField L Lbar` の座標の対を書こうとして

```lean
(fun q : (M : Type) × (M : Type) => (ψ q.1, ψ q.2))
```

と書くと、`(x : α) × β` が **`Sigma` の記法**なので
`Σ (M : Type), ↥M` にパースされ、`q.1 : Type`・`q.2 : q.1` になる。
エラーは

```
Application type mismatch: The argument q.fst has type Type of sort `Type 1`
but is expected to have type ↥M of sort `Type`
```

★対処: 型注釈の中では `↥M × ↥M` と書く(あるいは `M × M` で CoeSort に任せる)。
`(M : Type)` の形は `→+*` や `.baseChange` の引数では問題ない——`×` の左だけが危ない。

## `rw [← h]` で `h : addOrderOf X = l` を使うと motive が壊れる（2026-09-06）

`l` は点の型の中（`(range l).erase 0` の中の `veluQuotientFull` の点集合）にも
現れるので、`rw [← hQord]` は `l` を全部 `addOrderOf X` に置き換えようとして
`motive is not type correct` になる。
★対処: 前向きに書き換える。

```lean
have h := addOrderOf_nsmul_eq_zero X
rwa [hQord] at h        -- `addOrderOf X • X = 0` → `l • X = 0`
```


## `Skeleton` と `Found` の**同名の重複定義**を両方 `open` すると `Ambiguous term`(2026-09-06)

`ABC3.Skeleton.Divisor` と `ABC3.Found.Divisor` はどちらも
`IsCodimOnePt` / `PrimeDivisorPt` / `WeilDiv` / `IsNormalScheme` を持っている
(`decisions-pending.md` D8 の 3 番が解消を待っている重複)。両方 `open` すると

```
error: Ambiguous term
  IsNormalScheme
Possible interpretations:
  ABC3.Found.Divisor.IsNormalScheme X : Prop
  ABC3.Skeleton.Divisor.IsNormalScheme X : Prop
```

が名前の出現ごとに出る。★対処: **片方だけ `open`** し、もう片方は完全修飾で書く。
★なお 4 つとも **defeq** なので、`open` さえ分ければ
`Found` の定理に `Skeleton` の項をそのまま食わせられる(`exact` が通る。実測済み)。


## `ProfiniteGrp` の `M` に「点ごとの可換性」から `CommGroup` を付ける（2026-09-06、Skeleton/FrdI/Def28ProL で実測）

`Found/ProL/` は `variable {M : Type u} [CommGroup M] [TopologicalSpace M] …` で書かれているが、
`Skeleton` 側の statement は `(M : ProfiniteGrp.{u}) (_hcomm : ∀ a b : M, a * b = b * a)` である。
★`letI` で `Group` 構造に `mul_comm` を足すだけで橋が架かる（structure eta のおかげで
`CommGroup.toGroup` は元の `M.str` と defeq なので、`≃ₜ*` の型が変わらない）。

```lean
letI : CommGroup M := { (inferInstance : Group M) with mul_comm := _hcomm }
exact ⟨fun l => ABC3.Found.ProL.lPartGrp M l.1, …⟩
```

★もう一つの罠: `variable (M) in` が**付いていない**宣言（例 `isProL_lPartGrp`）では
`M` は implicit なので、`isProL_lPartGrp M l` と書くと

```
Application type mismatch: M has type ProfiniteGrp.{u} … but is expected to have type ℕ
```

（`M` が `l : ℕ` の位置に入る）。★対処は名前付き引数＋**coercion を明示**:
`isProL_lPartGrp (M := (M : Type u)) l.1`。`(M := M)` だけだと `CoeSort` が入らない。


## 触れない下流があるとき、statement に足す仮定は `∃` で述語の内側に畳む（2026-09-06、D8 で実測）

**失敗形**: `Skeleton` の `def` に仮定を足す（`ordAtDiv X x f` → `ordAtDiv X hnorm x f`）と、
その `def` を使う**別ファイルの述語**（`IsCartierDiv`）も引数が増え、
さらにその述語を使う**触ってはいけないファイル**（保留中の `Cartier/Theorem62.lean`）の
`IsCartierDiv X D` / `isCartierDiv_add X hD hE` が arity 不一致で落ちる。
★`[h : P]` の instance-implicit で逃げることはできない（`P` が class でないと
「invalid binder annotation」。class にしても下流で synthesize できない）。

**直し方**: 述語の側は引数を増やさず、仮定を **`∃` で内側に畳む**:

```lean
def IsCartierDiv (D : WeilDiv X) : Prop :=
  ∃ hnorm : IsNormalScheme X, ∀ x : X, …   -- 本文で hnorm を使う
```

★引数の数が変わらないので下流は無傷。しかも `hD : IsCartierDiv X D` から
`obtain ⟨hnorm, hD⟩` で仮定を**取り出せる**ので、`add` / `neg` 系の補題は
追加仮定なしで閉じる（`zero` 系だけが外から仮定を要る）。

★★**向きを間違えないこと**: `∀ hnorm, …` と書くと仮定が偽の枝が**空虚に真**になり、
「非正規スキームでは全部 Cartier（`cartierSubgroup = ⊤`）」という**新しい退化**を作る。
`∃` は逆に「仮定が偽なら述語も偽」なので**強める側**であり、退化を作らない。

★`hnorm : Prop` には definitional proof irrelevance が効くので、
2 つの `obtain` で別名の証明が出ても `ordAtDiv X hnorm₁ … = ordAtDiv X hnorm₂ …` は `rfl`。
`rw` が syntactic に噛み合わなくても最後の `exact` は通る。


## `↑(f⁻¹)`（`Units`）は `(↑f)⁻¹` と defeq ではない（2026-09-06）

**失敗形**: `show _ = ordPt X hnorm y ((f : K)⁻¹)` が

```
'show' tactic failed, pattern … is not definitionally equal to target
  -ordAtDiv X hnorm y ↑f = ordAtDiv X hnorm y ↑f⁻¹
```

で落ちる。`Units.mul` の `val` は `rfl` で割れる（`↑(f*g) = ↑f*↑g` は defeq）のに、
`Units.inv` の `val` は割れない（`f⁻¹.val` は構造体の `inv` 場であって `(f.val)⁻¹` ではない）。

**直し方**: `Units.val_inv_eq_inv_val` を先に当てる。`def` を挟んでいるときは
`simp only [ordAtDiv, Units.val_inv_eq_inv_val]` のように **`def` 名も `simp only` に入れて展開**
してから `rw [… ordPt_inv …]`。


## `Finsupp.single_eq_of_ne` の向きは直感と逆（2026-09-06）

**失敗形**: `single w 1 v = 0` を閉じようとして `Finsupp.single_eq_of_ne (Ne.symm hvw)` と書くと

```
Application type mismatch: hvw has type v ≠ w but is expected to have type w ≠ v
```

型は `(h : a' ≠ a) : (single a b) a' = 0` で、**添字が先ではなく後ろ**である。
さらに `simp [Finsupp.single_eq_of_ne hvw]` のように**部分適用**で渡すと
`b` がメタ変数のまま残り、linter に `This simp argument is unused` と言われる。

**直し方**: 素点の不等式そのものを渡すだけでよい（`simp [hvw]` / `simp [Ne.symm hvw]`）。
どちらが要るかは `single` の添字と評価点のどちらが左かで決まるので、
`simp` の unused 警告を見てから片方に決める。

## `nrRealPlaces` を `rw` すると `Fintype` の合成に落ちる（2026-09-06）

**失敗形**: `rw [← NumberField.InfinitePlace.card_real_embeddings]` の直後に

```
failed to synthesize instance of type class
  Fintype { φ // ComplexEmbedding.IsReal φ }
```

`IsReal` が decidable でないので `Subtype.fintype` が付かない。

**直し方**: 証明の冒頭に `classical` を置く。★同じ束で
`IntermediateField ℚ ℝ`（例: `ℚ⟮Real.sqrt 2⟯`）は `Real.instField` が noncomputable なので
`noncomputable abbrev` にしないと `dependsOnNoncomputable` で落ちる。
数体の実例を 1 つ作る手順は
`Check/FrdI/Thm64PicDegenerate.lean`（`sqrt2_isIntegral` → `NumberField` インスタンス →
`minpoly` → `finrank` → 実埋め込み → `card_add_two_mul_card_eq_rank` と `omega`）にある。

## 「同じ `j` の 2 曲線」と `ofJ` は mathlib に在る（2026-09-06、在庫）

`Mathlib/AlgebraicGeometry/EllipticCurve/ModelsWithJ.lean` と `.../IsomOfJ.lean`：

* `WeierstrassCurve.ofJ (j : F) [DecidableEq F] : WeierstrassCurve F`、`ofJ_j : (ofJ j).j = j`、
  無条件の `IsElliptic` インスタンスつき。`j ≠ 0, 1728` なら
  `ofJ_ne_0_ne_1728 : ofJ j = ofJNe0Or1728 j = ⟨j−1728, 0, 0, −36(j−1728)³, −(j−1728)⁵⟩`。
  `ofJNe0Or1728_c₄ = j(j−1728)³`・`ofJNe0Or1728_Δ = j²(j−1728)⁹`。
* `WeierstrassCurve.exists_variableChange_of_j_eq : E.j = E'.j → ∃ C, C • E = E'`
  ——★仮定は **`[IsSepClosed F]`**（分離閉体）である。数体の上では**そのままでは使えない**。
  有限次部分拡大へ降ろす配管（`C = ⟨u,r,s,t⟩` を含む `IntermediateField` を取る）が別に要る。

★★数体 `L` の素点 `p` で `v_p(j) < 0` なら、`ofJ j` は **拡大なしで** `p` で乗法還元になる
（`u = j` の変数変換で `v_p(c₄) = 0`）。証明は `Found/GenEll/VeluJExpNeg.lean`。

## 体の元には `linarith` が効かない（順序が無い）（2026-09-06）

`h : j - 1728 = 0` から `j = 1728` を出すのに `linarith` を使うと
`linarith failed to find a contradiction`。`L` は `LinearOrderedField` ではないからである。

**直し方**: `sub_eq_zero.1 h` / `linear_combination h` を使う。
`j + (-1728) = 0 → j - 1728 = 0` も `linear_combination h` で通る。

## `all_goals try ring` を保険で置くと linter が落ちる（2026-09-06）

`field_simp` が閉じてしまうと `'all_goals try ring' tactic does nothing` と
`this tactic is never executed` の 2 種の警告が出る。
**保険で置いた行は、1 回通ったら消す**（`sed -i '/all_goals try ring/d'`）。

## 大きな Lean ファイルを bash のヒアドキュメントで書くと落ちることがある（2026-09-06）

★数百行・非 ASCII 多数の内容を `cat > f.lean <<'EOF'` で書いたら
`/usr/bin/bash: -c: line 1: unexpected EOF while looking for matching \`''` で
**ファイルが 1 バイトも作られなかった**（同じ書き方で 60 行なら通る）。
**直し方**: 本体ファイルは Write ツールで書き、小さな差分だけ `sed -i` にする。

## `.absent` の本文に `--` を書くと `re:` の規約が消える（2026-09-06）

**失敗形**: `.absent "… 全体を grep、0 件(node tools/absent-recheck.mjs --try で再実行できる)。re:`Chebotarev|…`→0"`
と書いたのに、`node tools/absent-recheck.mjs` も `check.mjs` の G11 も
「パターンが無い」と数え続ける。

**原因**: 両者が使う `stripCommentsKeepStrings` は
`.replace(/--[^\n]*/g, …)` で **行コメントを文字列の中まで**潰す。
だから本文中の `--try` から行末までが空白になり、その後ろに書いた
`` re:`…` `` ごと消える。★`--brief` や `--json` を本文に書いても同じことが起きる。

**直し方**: `.absent` の本文に ASCII のハイフン 2 個を連ねない。
オプション名を挙げたいときは「`try` オプション」のように書く
（`——`（U+2014）や単独の `-`（日付の `2026-09-06` など）は安全）。

## `.absent` の `re:` にバックスラッシュを書くと、当たらない正規表現になる（2026-09-06）

**失敗形 1**: `` re:`AddMonoidHom\.toRealLinearMap` `` と書くと Lean が
`invalid escape sequence` で落ちる（`\.` は Lean の文字列エスケープに無い）。

**失敗形 2**（こちらが危ない）: 逃げて `` re:`AddMonoidHom\.toRealLinearMap` `` と書くと
Lean は通るが、`absent-recheck.mjs` の `readString` は
`\` を**そのまま 2 文字として**返す（unescape しない）ので、
取り出される正規表現は `AddMonoidHom\.toRealLinearMap`
＝「バックスラッシュ 1 個＋任意の 1 文字」になり、**永遠に 0 件**になる。
不在の主張が覆っても鳴らない、静かな穴である（実測で確認）。

**直し方**: `re:` の正規表現に ASCII のバックスラッシュを一切使わない。
`\.` は `[.]`、`\(` は `[(]`、`\s` は素の空白、`\b` は書かずに
`[A-Za-z0-9_.]*` などで前後を縛る。件数は
`node tools/absent-recheck.mjs --try '<正規表現>'` で書く前に測ること。

## `OpenSubgroup` を作るとき `rw [Subgroup.coe_map]` が当たらない（2026-09-06）

**失敗形**: `refine ⟨⟨H.map f, ?_⟩, ?_⟩` で `OpenSubgroup` を組むと、
第 1 の目標が `IsOpen (Subgroup.map f H).carrier` になる。ここに
`rw [Subgroup.coe_map]`（`↑(Subgroup.map f K) = f '' ↑K`）を当てると
"Did not find an occurrence of the pattern" で落ちる。構造体フィールドの
`isOpen'` は `SetLike` の `↑` ではなく `.carrier` で書かれているため。

**直し方**: `rw` の前に `show IsOpen ((… : Subgroup G) : Set G)` を挟んで
`SetLike` の coercion に戻す。`simp only [Subgroup.coe_map]` でも当たらない。

## `⁅a, b⁆` は `open scoped commutatorElement` が要る（2026-09-06）

**失敗形**: `have : x⁻¹ * y = ⁅a, b⁆ := by group` が
`failed to synthesize Bracket ↥H ↥H` で落ちる。記法自体は見えているのに
インスタンスが無い、という紛らわしい形で出る。

**直し方**: ファイル冒頭（`namespace` の中でよい）に
`open scoped commutatorElement`。なお `⁅a, b⁆ = a * b * a⁻¹ * b⁻¹` なので
`(h * y * h⁻¹)⁻¹ * y` は `⁅h, y⁻¹⁆` であって `⁅h⁻¹, y⁆` ではない
（`group` が閉じないときはここを疑う）。

## `Continuous (Subtype.val ∘ f.toFun)` に `rw [funext …]` は当たらない（2026-09-06）

**失敗形**: 部分群への写像の連続性を示すとき、`continuous_induced_rng.2` のあとの
目標は `Continuous (Subtype.val ∘ e.toFun)` の形で出る。
`have : (fun x => ↑(e x)) = fun x => … ; rw [this]` はパターンが合わずに落ちる。

**直し方**: `Continuous.congr` を使う——
`exact (連続性の証明).congr (fun x => (値の等式 x).symm)`。
`∘` と `.toFun` の差は defeq なので `congr` なら通る。

## `α (S.topologicalClosure) = (α S).topologicalClosure` は `Homeomorph.image_closure`（2026-09-06）

**失敗形**: `f.toHomeomorph.isOpenMap.image_closure_eq_closure_image` は**存在しない**
（`IsOpenMap` 側の名前を探しに行くと空振りする）。

**直し方**: `Subgroup` に落として `SetLike.coe_injective` してから
`f.toHomeomorph.image_closure _`（`⇑h '' closure s = closure (⇑h '' s)`）。
位相的アーベル化の移送（`Found/PGC/TopAbelianization.lean`）はこれ 1 行で足りる。

## `restrictNormalHom` と `restrictNormalHom_surjective` は `K₁`/`E` の役割が逆（2026-09-06）

**失敗形**: `AlgEquiv.restrictNormalHom` の暗黙引数は
`K₁ = 大きい体`・`E = 部分体`（`Gal(K₁/F) →* Gal(E/F)`）。ところが
`AlgEquiv.restrictNormalHom_surjective` のほうは
`K₁ = 部分体`・`E = 大きい体`（`(E : Type _)` が明示引数）。
同じ名前空間で規約が逆なので、片方に合わせて `(K₁ := …)` を書くと
`failed to synthesize Algebra 大 部分体` が出る。

**直し方**: 制限写像は `restrictNormalHom (F := F) (K₁ := Ω) (E := ↥E)`、
全射性は `restrictNormalHom_surjective (F := F) (K₁ := ↥E) Ω`。
`ker` を書くときは前者、全射性を渡すときは後者。

## `le_sup_left` は「項」なので包含に直接適用できない（2026-09-06）

**失敗形**: `h g (le_sup_left hg)` が
`Function expected at le_sup_left but this term has type ?a ≤ ?a ⊔ ?b` で落ちる。
`≤` を `∀ x ∈ …` に unfold して適用させたいのに、暗黙引数が決まらず
関数適用と読まれない。

**直し方**: 引数を明示する——`le_sup_left (a := H) (b := H') hg`。
（`Subgroup`・`IntermediateField` どちらでも同じ。）

## `jExp`／`.j` の中の曲線は `rw` で書き換えられない（2026-09-06、第 1445）

**失敗形**: `hZbc : Z.map f = veluQuotientFull (Y.map f) …` を
`rw [← hZbc] at h1M`（`h1M : jExp P' (veluQuotientFull (Y.map f) …) < 0`）で使うと
`motive is not type correct` になる。`jExp`（と `WeierstrassCurve.j`）は
`[W.IsElliptic]` を**曲線に依存するインスタンス引数**として取るからである。

**直し方**: `rw` を使わず `jExp_congr_j P' _ _ (j_congr_curve hZbc.symm)` で
**`jExp` の値どうしの等式**を作ってから、それを `rw` する（ℤ の等式なので安全）。
☆`j_congr_curve : W = X → W.j = X.j`（`Found/GenEll/JScale.lean`）が
`congrArg WeierstrassCurve.j` の代用である（`congrArg` は依存型なので通らない）。
★同じ理由で `E.baseChange A` と `E.map (algebraMap _ A)` の往復も `rw` ではなく
`have h' : … := h`（defeq に任せる）で行うこと。

## `by_contra` のあと `push_neg` を `not_le.mp` に置き換えたら whnf timeout（2026-09-06、第 1445）

**失敗形**: `by_contra hcon; push_neg at hcon` は通るのに、
`by_contra hcon0; have hcon : jExp p E' < 0 := not_le.mp hcon0` に書き換えると
**証明全体が `(deterministic) timeout at whnf`** で落ちた（`set` で置いた
局所定義 `E'` を展開しに行くらしい）。

**直し方**: `push_neg` のままにする（deprecation 警告は無害）。
どうしても消したいなら `set_option maxHeartbeats` を上げる。

## `restrictNormalHom` を**主張の中**で使うには `Normal` が instance でなければならない（2026-09-06、Λ3 位相版）

**失敗形**: `normal_lubinTateClosure`（在庫の**定理**、explicit 引数 `K hq hπmax …`
を取る）を `haveI` で入れても、`theorem foo … : … restrictNormalHom … = …` の
**statement の elaboration は `haveI` より前**に走るので
`failed to synthesize Normal K.carrier ↥(lubinTateClosure …)` で落ちる。
explicit 引数を持つ定理は `attribute [local instance]` にもできない。

**直し方**: そのセクションで
`variable [Normal K.carrier (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf)]`
と**インスタンス変数として受け取る**。`Normal` は `Prop` クラスなので
どのインスタンスを渡しても命題は同じ。仮説なしの形（最終定理）で
`haveI := normal_lubinTateClosure …` を置けば消える。
☆同じ形は `IsGalois`・`FiniteDimensional` を要求する statement でも起きる。

## 二重の第一同型定理で作った `MulEquiv` の値は `simp only [定義名, MulEquiv.trans_apply]` でほどける（2026-09-06、Λ3 位相版）

**失敗形と思い込んでいたもの**: `noncomputable def Φ := by haveI …; exact (e₁.symm.trans (e₂.trans …))`
の形で作った同型は `haveI`（= `let_fun`）が挟まるので値の計算ができない、と考えて
別定義を作り直しそうになる。

**実際**: `simp only [Φ, MulEquiv.trans_apply]` で素直に
`e₄ (e₃ (e₂ (e₁.symm x)))` までほどける。あとは
`(QuotientGroup.quotientKerEquivOfSurjective φ hφ).symm (φ g) = QuotientGroup.mk g`
（`MulEquiv.symm_apply_eq` → `rfl`）と `QuotientGroup.quotientMulEquivOfEq_mk` を
`simp` に渡し、最後は `rfl`。

## `↑(GrpCat.of (Multiplicative ℤ))` の元に `Multiplicative.toAdd` を当てられない（2026-09-06、Λ5 `Ẑ`）

**失敗形**: `FiniteIndexNormalSubgroup ↑(GrpCat.of (Multiplicative ℤ))` の要素 `a` に
`Multiplicative.toAdd a` と書くと
`failed to synthesize HPow … ↑(GrpCat.of (Multiplicative ℤ)) …`。
型注釈 `(a : Multiplicative ℤ)` も**効かない**（Lean が coe を探しに行く）。

**直し方**: `show T from a` で defeq を強制する。

```lean
obtain ⟨n, rfl⟩ : ∃ n : ℤ, (Multiplicative.ofAdd n : Multiplicative ℤ) = a :=
  ⟨Multiplicative.toAdd (show Multiplicative ℤ from a), rfl⟩
```

以後 `a` が `Multiplicative.ofAdd n` の形になるので、`^`・`toAdd` の
インスタンスが素の `Multiplicative ℤ` 側で解決する。
☆同じ理由で `rwa [heq] at hpow` が「パターンが見つからない」と言うことがある
（`^` が圏の対象の coe 側で elaborate されている）。`exact heq ▸ hpow` は通る。

## `dif_neg h` は coe の下の `dite` に届かないことがある（2026-09-06、Λ5 `Ẑ`）

**失敗形**: `def f (N : ℕ) := if h : N = 0 then 0 else (choose …)` の仕様を
`simp only [f, dif_neg hN]` で開こうとすると、`↑(if h : N = 0 then …)` の
ように **coe の下に残った `dite`** だけ書き換わらず型不一致になる。

**直し方**: 先に `theorem f_of_ne (hN : N ≠ 0) : f N = (choose …) := dif_neg hN` を
作り、`rw [f_of_ne K hN]` で当てる。`rw` は defeq 込みで一斉に置き換わる。

## `ProfiniteCompletion.lift` の成分計算は `IsLimit.fac` を掘らずに `lift_unique` で出る（2026-09-06、Λ5 `Ẑ`）

**やりたいこと**: `ψ := ProfiniteGrp.ProfiniteCompletion.lift f : completion G ⟶ P` の
開正規部分群 `U` での成分 `proj U (ψ x)` を `x` の成分で書く。
`lift` は `P.isLimitCone.lift ⟨…⟩` と**その場で作った錐**で定義されているので
`IsLimit.fac` を使おうとすると錐を再現する羽目になる。

**効いた直し方**: 両辺とも `completion G ⟶ (diagram P).obj U` なので
`lift_unique`（`eta` を前置して一致すれば等しい）で 3 行。

```lean
theorem lift_comp_proj (f : G ⟶ GrpCat.of P) (U : OpenNormalSubgroup P) :
    lift f ≫ ProfiniteGrp.proj U
      = (ProfiniteGrp.limitCone (diagram G)).π.app (preimage f U)
        ≫ ProfiniteGrp.ofFiniteGrpHom (quotientMap f U) := by
  apply lift_unique
  rw [Functor.map_comp, ← Category.assoc, lift_eta]
  ext y
  rfl
```

あとは `ConcreteCategory.congr_hom … x` と `simp only [ProfiniteGrp.comp_apply]`。
`quotientMap f U` は（`preimage f U = f⁻¹ U` なので）**常に単射**なので、
`ψ x = 1` から `x.val (preimage f U) = 1` が出る。
`ProfiniteGrp.limit` の元は `Subtype` なので `x.val j` / `x.2 h.hom`（整合性）が使える。

## `Skeleton` の主張に `Found` の証明を配線すると循環する ── 定義を 3 層目に降ろす（2026-09-06、pGC Prop 1.2 / Cor 1.3）

**症状**: `Skeleton/X.lean` に「対象の定義 + 主張」が同居していて、`Found` の証明ファイルが
その**定義**を使っている。主張の `sorry` を `Found` の定理で埋めようとすると
`Skeleton の主張 → Found の証明 → Skeleton の定義` で import が循環する。
（実例: `Skeleton/PGC/Section1Cor13.lean` の `inertia` を `Found/PGC/UnramifiedCriterion.lean`
以下の PGC の Found スタック 66 本が使っていた。）

**直し方**: 定義だけを新しい `Skeleton/…Defs.lean` に降ろして 3 層にする。

```
Skeleton/PGC/Section1Defs.lean   （定義。import は Setup + Interface だけ）
   ↑                                   ↑
Found/PGC/…（証明）              Skeleton/PGC/Section1.lean（主張。Found へ委譲）
```

`Found` 側の import を新ファイルへ付け替えるだけで済む（1 行）。
定義の中身は 1 文字も変えない。★循環の有無は木全体で 1 回測ること
（`import` を辿って DFS。2,167 モジュールで 1 秒）。

**ついでに踏む罠**: `check.mjs` の **G8** は「`sorry` 0 の Skeleton 定理」に対し、
**`Found/` 側に同じ `[paper] item` の `.src` が在ること**を要求する。
`Skeleton` で証明を書き切ると G8 が鳴って NG が増える。
`Found` に無条件版の定理と `def <名前>.src : ABC3.Meta.Source := { … item := "Proposition 1.2" … }`
を置き、`Skeleton` は 1 行で委譲する（`Skeleton/GenEll/Section4.lean` の Lemma 4.1 と同じ形）。
`.needs` は Skeleton の theorem/lemma にしか要求されないので `Found` 側には要らない。

## `namespace ABC3.Found.NF` の中の `open NumberField` は `ABC3.Found.NumberField` を開く（2026-09-06、判断 D11）

`Found/NumberField/RatHeightOne.lean` を書き始めたとき、`namespace ABC3.Found.NF` の中で

    open NumberField IsDedekindDomain

と書いたら `𝓞` が **Unknown identifier** になった（`autoImplicit=false` なので
「暗黙束縛にもできない」と言われる）。★原因は名前解決である ——
`ABC3.Found.NumberField` という名前空間が実在する（`Found/NumberField/Theory.lean` と
`PrimeDivisorsOfValues.lean` が使っている）ので、`ABC3.Found.NF` の中の `NumberField` は
**そちらに解決され**、mathlib の `NumberField` は開かれない。

★直し方は `_root_.` を付けるだけ。既存の `Found/NumberField/*.lean` は全部そう書いてある:

    open _root_.NumberField IsDedekindDomain Ideal
    open scoped _root_.NumberField

★`open scoped` の側にも `_root_.` が要る（`𝓞` は scoped notation なので、
片方だけ直すと記号だけ出てこない）。

## `def f (_x : A) : B := sorry` は**定数関数**に展開される（2026-09-06、判断 D16）

`Skeleton/Divisor/ArithDivisor/Example63.lean` の

    noncomputable def degArith (_x : ArithPhiGp L) : ℝ := sorry

について、`degArith L x = degArith L y` が **`rfl` で通る**。
`sorry` の項は引数に依存しないからである。したがって

    theorem foo : Function.Surjective (degArith L)

は現在の環境で**反証可能**であり（`0 = 1` が出る）、この `sorry` は
`degArith` に本体を与えない限り**原理的に閉じない**。

★これは place-holder についての主張であって原典についての主張ではない
（反証は `degArith` 経由で `sorryAx` に依存する）。
★★**使い道**: 「`sorry` 本体の `def` の非定数性（全射・単射・非退化）を主張する
skeleton」は、その場で機械的に反証できる。退化の検査を書く前に 3 行で判る。

## `rw` で曲線を置き換えると `jExp` の motive が壊れる（2026-09-06、第 1446）

`jExp p W` は `[W.IsElliptic]` を**引数に取る**ので、

    have heq : C • X = Y
    rw [← heq] at hres          -- hres : jExp p Y < 0

は `motive is not type correct` で落ちる（`fun _a => jExp p _a < 0` を作れない：
`_a.IsElliptic` のインスタンスが `_a` に依存する）。

★直し方は**等式を `jExp` の外に出す**——`j_congr_curve`（曲線が等しければ `j` が等しい）と
`jExp_congr_j`（`j` が等しければ `jExp` が等しい）を並べて `omega` で閉じる:

    have h1 : jExp p (C • X) = jExp p Y := jExp_congr_j p _ _ (j_congr_curve heq)
    have h2 : jExp p (C • X) = jExp p X := jExp_variableChange p _ C
    omega

☆同じ理由で `minDeltaExp`・`SemistableAt` を含む書き換えも詰まる。
`SemistableAt` はインスタンスを取らないので `rw` が通る——**`jExp` だけが違う**。

## `jExp` を結論に持つ定理は `IsElliptic` を**仮引数に置く**（2026-09-06、第 1446）

    theorem foo ... : jExp p (veluQuotientFull E S) < 0

は statement の elaboration で `(veluQuotientFull E S).IsElliptic` を探しに行き、
`synthInstance failed` で落ちる（証明の中で `haveI` しても手遅れ）。
★`[hVell : (veluQuotientFull E S).IsElliptic]` を**仮引数に足す**。
☆`SemistableAt p (veluQuotientFull E S)` を結論にする版は
インスタンスを要らないので、`jExp` 版へ移すときにこれが必ず出る。

## `↥(myDef)` に mathlib の `↥v.integer` 補題を `rw` できない（2026-09-06、Λ5）

    noncomputable def OK : ValuationSubring F := (Valued.v : Valuation F NNReal).valuationSubring
    example (w : ↥(OK)) : ¬ IsUnit w ↔ ‖(w : F)‖ < 1 := by
      rw [Valuation.Integer.not_isUnit_iff_valuation_lt_one]   -- ✗ pattern not found

`rw` は instance-reducible までしか展開しないので、`↥(OK)` は `↥v.integer` と
**構文的に**一致せず、`{x : ↥v.integer}` の暗黙引数が埋まらない。

★直し方は**補題を項として適用し、引数で型を合わせる**（既定透明度なので `def` が展開される）:

    have h := Valuation.Integer.not_isUnit_iff_valuation_lt_one
      (v := (Valued.v : Valuation F NNReal)) (x := w)
    rw [valuation_eq_nnnorm] at h
    exact h.trans ⟨fun hh => by exact_mod_cast hh, fun hh => by exact_mod_cast hh⟩

## 完備化に `Valued` が二重に付く（2026-09-06、Λ5）

底 `F` に `instance : Valued F NNReal := NormedField.toValued` を置くと、
mathlib の `Valued.valuedCompletion`（priority 既定の instance）が
`UniformSpace.Completion F` に**自動で** `Valued` を作る。そこへ自分でも
`NormedField.toValued` を置くと菱形になり、`Valued.v z = ‖z‖₊` が `rfl` で通らなくなる。

★片側に寄せる。底では `Valued` を **instance にせず**（必要な
`CompletableTopField` は `letI : Valued F NNReal := NormedField.toValued; Valued.completable`
で借りる）、完備化側だけに `NormedField.toValued` を置く。こうすると
`Valued.v z = ‖z‖₊` が `rfl` になり、以降の付値の議論がすべてノルムの議論に落ちる。

## `f z` の値を 2 層剥がす `rfl` は kernel を止める ── `_apply` 補題を挟む（2026-09-06、Λ5b）

    noncomputable def F : A →+* ↥(myValuationSubring) where
      toFun z := ⟨⟨(z : Base), h₁⟩, h₂⟩
      ...
    theorem F_coe (z : A) : ((F z : ↥S) : Mid) = ... := rfl   -- ✗ (kernel) deterministic timeout 57 秒

`RingHom` の構造体の展開と `Subtype.val` を **2 層**同時に要求すると kernel が落ちる
（#59 と同型。elaborator は 0.1 秒で通るのに kernel だけが止まるので気付きにくい）。
`w : ↥S` を変数にした `(w : Mid) = w.1 := rfl` は 0.03 秒で通る ── 詰まるのは
**構造体インスタンスの適用と Subtype の反復**の組み合わせである。

★直し方は**トップレベルの `rfl` で `_apply` 補題を先に置き、あとは `rw` で 1 層ずつ**:

    theorem F_apply (z : A) : F z = ⟨⟨(z : Base), h₁⟩, h₂⟩ := rfl     -- 0.06 秒
    theorem F_coe (z : A) : ((F z : ↥S) : Mid) = ⟨(z : Base), h₁⟩ := by rw [F_apply]
    theorem norm_F (z : A) : ‖((F z : ↥S) : Mid)‖ = ‖z‖ := by rw [F_apply]; rfl

`_apply` は両辺が同じ項なので kernel の仕事が消え、`rw` のあとに残る `rfl` は 1 層だけになる。

## `set x := e with h` は他の仮説の**型**まで書き換えて別の変数を生む（2026-09-06、Λ5b）

    (σ : ↥K⟮x⟯ ≃ₐ[F] ↥K⟮x⟯) (hσ : P σ)
    set L := K⟮x⟯ with hL      -- ここで σ : Gal(↥L/F) と σ✝ : Gal(↥K⟮x⟯/F) に**分裂する**
    rw [hσ]                    -- ✗ pattern not found（hσ は σ✝ の話）

`set` は目標だけでなく文脈の型にも置換をかけるため、束縛済みの変数が
「置換後の型を持つ新しい変数」として複製され、古い仮説は `σ✝` を指したままになる。

★長い型を短くしたいだけなら `set` を使わず**そのまま書く**か、`set` は
仮説に現れない局所的な項（構成した `⟨_, _⟩` など）だけに使う。

## `rw [← map_add]` は「もう片方の辺が `def` のまま」だと当たらない（2026-09-06、D10 の配線）

**失敗形**: `arithDivOfElt L (f * g) = arithDivOfElt L f + arithDivOfElt L g` を
`show (equiv.symm _) = _` で左辺だけ開いてから `rw [← map_add]` すると

```
Did not find an occurrence of the pattern ?f ?x + ?f ?y
in the target expression
  equiv.symm ⟨arithDiv (f * g), ⋯⟩ = arithDivOfElt L f + arithDivOfElt L g
```

`show` は**書いた側しか**開かないので、右辺は `def` の名前のまま残り
`?f ?x + ?f ?y` に一致しない。

**直し方**: `unfold <def名>` で**両辺**を開いてから `rw [← map_add]`。
`simp only [<def名>]` でもよいが、`simp only` は続けて余計な正規化をかけて
`← map_add` の形を壊すことがあるので `unfold` の方が安全である。

## `Finsupp.single_apply` は素点の上では `classical` が要る（2026-09-06）

**失敗形**: `rw [Finsupp.single_apply]` が

```
failed to synthesize instance of type class Decidable (v = w)
```

`InfinitePlace L` / `FinitePlace L` / `Sum` に `DecidableEq` が付いていないためで、
補題が無いのではない。**定理の先頭に `classical` を置くだけ**で通る
（`by_cases` + `Finsupp.single_eq_of_ne` に逃げると今度は向き（上の項）で嵌まる）。


## ★文字列の置換をヒアドキュメントの Python で繰り返さない(2026-09-06、本体セッションが 6 往復浪費)

**失敗形**: `tools/source-text.py` の 2 行を直すのに、
`python - <<'PYEOF'` の中で `s.replace(old, new)` を書いて `assert s.count(old) == 1` で
**6 回連続で外した**。原因は毎回違った:
インデントの実体差 / `\n` の層の数え違い / 1 つ目の assert で止まって 2 つ目が未適用 /
行番号の 0-based と 1-based。

**直し方**: ★★**Read してから Edit ツールを使う**。一発で通った。
Edit は実ファイルの文字列をそのまま受けるので、
**シェル→Python →正規表現の 3 層のエスケープが消える**。

★**往復回数がそのまま費用**である(2026-09-06 の実測: 費用の 99.36% は cache_read、
1 往復 547K token)。**2 回外したら手法を変える**こと。

★例外: 新規ファイルの作成や、同じ置換を多数のファイルに当てるときは Python が正しい。
問題は「1 ファイルの数行を直す」のに Python を使うことである。


## `ProfiniteGrp.limit` の成分は `G ⧸ H` と defeq だが syntactic には別（2026-09-06、経路 Λ9）

**失敗形**: `x : ZHat`（= `ProfiniteCompletion.completion (GrpCat.of (Multiplicative ℤ))`）の
成分 `x.val H` の型は `((diagram G).obj H).toProfinite.toTop` で、`G ⧸ H.toSubgroup` と
**defeq だが syntactic には別**。そのため

* `rw [QuotientGroup.mk_pow]` / `rw [← limit_pow_val]` が
  `Did not find an occurrence of the pattern` で落ちる
  （末尾に `The target expression is not type-correct under the instances transparency level` が付く）
* `exact mul_comm _ _` が `CommMagma ↑((diagram G).obj H)` を探して落ちる

**直し方**: 3 つとも「型を書く」で通る（0.05 秒）。

1. **`have` の型に書く**: `have h : ((x.val H : _) : Multiplicative ℤ ⧸ H.toSubgroup) ^ m = 1 := h1`
   —— `have`/`show` の型注釈は defeq で通り、以後 `rw` が効くようになる。
2. **`obtain` の型に書く**: `obtain ⟨a, ha⟩ : ∃ a : Multiplicative ℤ, (QuotientGroup.mk a : _ ⧸ _) = ... := QuotientGroup.mk_surjective _`
   —— `a : ↑(GrpCat.of (Multiplicative ℤ))` のまま取ると
   `Multiplicative.toAdd a` が `Multiplicative ↑(GrpCat.of ...)` で誤って通り、
   後の `toAdd_pow` が全部落ちる。
3. **★項の型注釈 `(e : T)` は「落ちる」**。`exact mul_comm (x.val H : Multiplicative ℤ ⧸ H.toSubgroup) _`
   は依然 ProfiniteGrp 側の型で探しにいく。**名前つき暗黙引数で型を固定する**こと:
   `exact @mul_comm (Multiplicative ℤ ⧸ H.toSubgroup) _ (x.val H) (y.val H)`。

★`ProfiniteCompletion.completion` は `abbrev` ではなく `def` なので、
`limit_pow_val (diagram G) x H k` のような**補題の適用（`exact`）は通る**（default transparency）が、
`rw` は通らない。**補題は `ZHat` 版に言い換えて置く**（`zhat_pow_val`）のが安い。

## #70 `div_le_iff₀` / `div_lt_iff₀` が作る積の**左右**を当てにいかない（2026-09-06、Λ6-M1）

**失敗形**: `rw [div_le_iff₀ h] at hkey` で `‖z‖ ≤ ‖π‖ ^ m * ‖π‖` になるか
`‖π‖ * ‖π‖ ^ m` になるかは、周りの `have` の書き方ひとつで入れ替わる。
補助等式 `hexp : ‖π‖ ^ m * ‖π‖ = ‖π‖ ^ (m+1)` を `rw [hexp] at hkey` で当てると
半分の確率で `Did not find an occurrence of the pattern` になり、
`mul_comm` を足すか外すかで**同じ往復を 2 回**やる羽目になる。

**直し方**: 積の向きを当てにいかず、**等式を linarith/nlinarith に渡す**
（`linarith [hexp, hkey]`）か、`rw [hexp]` の代わりに
`calc`／`linear_combination` を使う。どうしても `rw` にするなら
`mul_comm` を**一度だけ**入れて、外れたら向きを変えるのではなく手法を変えること。

## #71 `hcard ▸ orderOf_dvd_natCard a` は motive が壊れる（2026-09-06、Λ5b′）

**失敗形**: `Nat.card (Multiplicative (ZMod N)) = N` を `hcardG` として
`have hdvd : orderOf a ∣ N := hcardG ▸ orderOf_dvd_natCard a` と書くと、
`▸` が **`a` の型の中の `ZMod N`** まで一緒に書き換えようとして
`orderOf (a : Multiplicative (ZMod (Nat.card ...)))` という化物を作る
（`rw [← hcardG]` にしても `motive is not type correct`）。
`Nat.card G` の `G` が引数の型に現れる補題（`orderOf_dvd_natCard`,
`Nat.card_zpowers`, `Subgroup.eq_top_of_card_eq` 系）で必ず踏む。

**直し方**: ゴール側を触らず、**補題を仮説に落としてから仮説を書き換える**:

```lean
have h := orderOf_dvd_natCard a
rwa [hcardG] at h    -- 型の中の ZMod N は動かない
```

## #72 `map_mul` を `rw` で数え間違えると `whnf` timeout に化ける（2026-09-06、Λ7 Lemma 4.11）

**失敗形**: 準同型 `φ` について `φ (a * b * a⁻¹ * b⁻¹) = 1` を出そうとして
`rw [map_mul, map_mul, map_inv, map_inv]` と書いた。積は **3 個**なので
`map_mul` が 1 回足りず、ゴールは `φ (a * b) * (φ a)⁻¹ * (φ b)⁻¹ = 1` のまま。
次の `exact` が `φ (a*b)` と `φ a * φ b` を合わせにいって
`(deterministic) timeout at 'whnf', maximum number of heartbeats (200000)`。
★エラーが「rw が外れた」ではなく **timeout** として出るので、原因が
`rw` の回数だと気づくまでに往復を無駄にする（実測 11 秒 × 1）。

**直し方**: `rw` で `map_mul` を数えない。**`simp only [map_mul, map_inv]` を
`have h : φ (…) = 1 := by …` の中で使い**、外側は `exact h` にする。

## #73 `x ∈ H ⊓ K`（Subgroup）を `refine ⟨?_, ?_⟩` で割ると `toSubmonoid` が露出する（2026-09-06、Λ7）

**失敗形**: `⊢ g ^ m ∈ L.fixingSubgroup ⊓ Km.fixingSubgroup` に
`refine ⟨?_, hpow g⟩` を当てると、第 1 のゴールが
`g ^ m ∈ ↑(AlgEquiv.restrictNormalHom ↥L).ker.toSubmonoid` という
**部分モノイドの coe** になり、続く `rw [MonoidHom.mem_ker]` が
`Did not find an occurrence of the pattern ?x ∈ MonoidHom.ker ?f` で落ちる。

**直し方**: `Subgroup.mem_inf.mpr ⟨h1, h2⟩` を使う（`h1`・`h2` は
`x ∈ H`・`x ∈ K` を**別の `have` で作っておく**）。`⊓` を anonymous
constructor で割らない。

## #74 `IsHausdorff.haus inferInstance` はイデアルが決まらず stuck（2026-09-06、Λ6 Dwork 加法版）

**失敗形**: `refine IsHausdorff.haus inferInstance _ (fun N => ?_)` が

    typeclass instance problem is stuck, it is often due to metavariables
      IsHausdorff ?m.243 ↥(unramifiedCompletionInt K)

`IsHausdorff.haus` の結論は `x = 0` で、イデアル `I` は**結論に現れない**ので
暗黙引数が埋まらない。同じ `IsAdicComplete` から取る `IsPrecomplete.prec` は
Cauchy 条件 `hcauchy` の型に `I^m • ⊤` が出るので**そちらは埋まる**——
片方だけ落ちるので原因に気づきにくい。

**直し方**: `(I := …)` を明示する。

    refine IsHausdorff.haus (I := IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K))
      inferInstance _ (fun N => ?_)

★★**この罠は 2 つ目の顔を持つ**（2026-09-08、pGC §2 分岐入力）。`(I := …)` を書いても
`inferInstance` を**書き忘れる**と、今度は型不一致になる（`haus` は構造体フィールドなので
`self` を明示的に取る＝引数が 1 つ多い）:

```
error: Type mismatch
  IsHausdorff.haus ?m.109 ?m.110
has type
  (∀ (n : ℕ), ?m.110 ≡ 0 [SMOD maximalIdeal ↥𝒪[K.carrier] ^ n • ⊤]) → ?m.110 = 0
but is expected to have type
  ↑u - 1 = 0
```

★**引数は `(I := …) inferInstance x h` の 4 つ**である。`(M := …)` を足しても

```
error: typeclass instance problem is stuck
  Module ↥𝒪[K.carrier] ?m.101
```

が出るだけで直らない。

## #75 `SModEq` は `I^n • ⊤` の形なので、差の所属に落とす補題を 1 本置く（2026-09-06、Λ6）

`IsPrecomplete.prec` / `IsHausdorff.haus` が使う形は `f m ≡ f n [SMOD I ^ m • ⊤]` で、
環を自分自身の加群と見た `I ^ m • (⊤ : Submodule R R)` である。往復のたびに
`SModEq.sub_mem` と `I • ⊤ = I` を手で当てると読めなくなるので、1 本だけ置く:

    theorem sModEq_pow_iff {R : Type*} [CommRing R] (I : Ideal R) (n : ℕ) (x y : R) :
        x ≡ y [SMOD I ^ n • (⊤ : Submodule R R)] ↔ x - y ∈ I ^ n := by
      rw [SModEq.sub_mem]; simp

★`I • (⊤ : Submodule R R) = I` は `simp` が一発で閉じる（専用の補題名を探さなくてよい）。

## #76 `choose` で漸化式を組むなら存在文は `∃ t, (仮定 → 結論)` にする（2026-09-06、Λ6 Dwork 乗法版）

**失敗形**: 1 段補題を素直に

    theorem step (n : ℕ) (w : Oˣ) (hw : ↑w - 1 ∈ (π^n)) : ∃ t : Oˣ, P n w t

と書いて `choose T hT using step` すると、`T` の型が

    T : ∀ (n : ℕ) (w : Oˣ), ↑w - 1 ∈ (π^n) → Oˣ

になり、**証明を引数に取る関数**なので `Nat.rec` で数列 `w_{n+1} := f (T n (w n)) (w n)`
が組めない（`w n` が仮定を満たす証明は帰納法の途中でしか手に入らない）。

**直し方**: 仮定を存在の**内側の含意**に押し込んで、結論だけを条件つきにする。

    theorem step (n : ℕ) (w : Oˣ) : ∃ t : Oˣ, (↑w - 1 ∈ (π^n) → P n w t)

こうすると `T : ℕ → Oˣ → Oˣ` と `hT : ∀ n w, ↑w - 1 ∈ (π^n) → P n w (T n w)` に割れる。
仮定が `t` に依存しないので中身は変わらない。証明側は `by_cases` で
「仮定が成り立つなら本来の構成、成り立たないなら `1`」と分ければよい。
★段の場合分け（`n = 0` は剰余体の `q-1` 乗根、`n ≥ 1` は加法版）もこの形の中に隠せる。

## #77 `Function.Surjective (fun x => f x)` の目標に `rw` は当たらない（2026-09-06、Λ6）

**失敗形**: `intro u; refine ⟨ξ, ?_⟩` のあと目標が

    (fun ξ => (unramGalCompletionUnits K σ) ξ * ξ⁻¹) ξ = u

という **β-簡約されていない形**のままなので、`rw [hξ]` が
「`(unramGalCompletionUnits K σ) ξ` が見つからない」で落ちる。

**直し方**: `show` で β-簡約した形を書いてから `rw` する（`beta_reduce` / `simp only []`
でもよいが `show` が一番安い）。

    show unramGalCompletionUnits K σ ξ * ξ⁻¹ = u
    rw [hξ, ...]

## #78 `simpa using h` が linter に叱られたら `sub_sub_cancel` を疑う（2026-09-06、Λ6）

`x ∈ 𝔪` と `x - y ∈ 𝔪` から `y ∈ 𝔪` を出すとき `Submodule.sub_mem _ hx hxy` は
`x - (x - y) ∈ 𝔪` を返す。`simpa using` で通るが
`Try 'simp at h' instead`（= simp が目標まで閉じている）という警告が出て意図が濁る。
`rwa [sub_sub_cancel] at h` と書けば `a - (a - b) = b` の 1 手だと読める。

## #79 「`N` を法として等しい 2 元」の乗り換えは商群を 1 回作るのが最短（2026-09-06、N13）

**場面**: Frobenius の持ち上げ `σ₀` を「`F` を固定する `σ`」に取り替えたとき、
`σ₀` に課された条件（`σ₀^k ∈ N ↔ m ∣ k` など）を `σ` へ移したい。
`σ₀⁻¹σ ∈ N`（`N` 正規）だけが手掛かり。

**失敗形**: `(σ₀ n)^k σ₀^{-k} ∈ N` を `ℤ` の帰納法で直接示そうとすると長い。

**直し方**: `QuotientGroup.mk' N` を通す。3 行で済む。

    have hq : (QuotientGroup.mk' N) σ = (QuotientGroup.mk' N) h := by
      simp only [QuotientGroup.mk'_apply]; exact QuotientGroup.eq.mpr hmem
    -- あとは (QuotientGroup.eq_one_iff _) と map_zpow で往復するだけ

★`QuotientGroup.mk'_eq_mk'` は**無い**（`Unknown constant` になる）。
`QuotientGroup.eq : (↑a = ↑b) ↔ a⁻¹ * b ∈ s` を `mk'_apply` で剥がしてから使う。
★`⟨σ⟩ ⊔ N = ⊤` の乗り換えだけは商群も正規性も要らない
（`σ = h * (σ⁻¹h)⁻¹ ∈ ⟨h⟩ ⊔ N` で `Subgroup.zpowers_le` に流す）。
★ついでに「`_root_.mem_fixingSubgroup_iff` は M が明示引数」の罠は `IntermediateField.mem_fixingSubgroup_iff` にもある
（`K` と `σ` が**両方とも明示引数**：`(IntermediateField.mem_fixingSubgroup_iff F g).mp h`。
`_` を 1 つにすると `Invalid field 'mp' ... Function.mp` という無関係なエラーが出る）。

★**逐語のエラー文**（2026-09-08、pGC 中心化群。★この節は 2026-09-07 まで逐語を持たず、
同じ罠に別の顔で 2 度落ちている）。第 1 引数（中間体）を省いて `g` を先に渡すとこうなる:

```
error: Application type mismatch: The argument
  g
has type
  Gal(E/F)
but is expected to have type
  IntermediateField ?m.53 ?m.55
in the application
  IntermediateField.mem_fixingSubgroup_iff g
```

```
×  (IntermediateField.mem_fixingSubgroup_iff g).mp hg
○  (IntermediateField.mem_fixingSubgroup_iff L g).mp hg y hy
○  rw [IntermediateField.mem_fixingSubgroup_iff] at hg   -- rw なら引数を書かなくてよい
```

★`rw ... at` 版は引数順を知らなくても通るので、**まず rw で開く**のが安い。


## #80 `lean_check` の断片で `open X in section … end` と書くと open が中身に届かない（2026-09-06、Λ6b）

**症状**: MCP の `lean_check` に

    open ABC3.Found.PGC in
    open scoped Valued in
    section
    variable ...
    #check (𝒪[K.carrier] : Type)
    end

を投げると `Unknown identifier `𝒪`` が出る。`𝒪[…]` は `open scoped Valued` の
スコープ記法なので「開いたはず」なのに届いていない。

**原因**: `open ... in` は**直後の 1 コマンドだけ**を修飾する。`section` はそれ自体が
1 つのコマンドなので、`in` は `section` に掛かって終わり、中身には届かない。

**直し方**: `in` を外して、`section` の**中**に `open` を並べる。

    section
    open ABC3.Found.PGC
    open scoped NNReal Valued
    variable {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    ...
    end

★`#check @Foo.bar` のように完全修飾で書いた行だけは通ってしまうので、
「一部だけ通って一部が Unknown identifier」という**紛らわしい**出方をする。

## #81 REPL で `open` を忘れると「Unknown identifier」ではなく **65 秒の heartbeat 焼き**になる（2026-09-06、Λ6 DworkTheta）

**症状**: `lean_check` に `Found/PGC` 用の断片を

    namespace ABC3.Found.PGC
    variable {p : ℕ} [Fact p.Prime]
    noncomputable def foo (K : PAdicLocalField p) …

と投げると、`Function expected at PAdicLocalField` に続いて
`(deterministic) timeout at isDefEq, maximum number of heartbeats (1000000)`
が出て **65 秒**返ってこない。`set_option maxHeartbeats 1000000` を付けていると
その分だけ焼く。

**原因**: `PAdicLocalField` は `ABC3.Skeleton.PGC` にある。REPL は
`autoImplicit` が有効（#62）なので、未知の識別子は**型未定のメタ変数**として
自動束縛され、その先の `rfl` / `isDefEq` がメタ変数を相手に探索し続ける。
「Unknown identifier」で即座に落ちてくれない。

**直し方**: `Found/PGC` の断片には必ず

    namespace ABC3.Found.PGC
    open ABC3.Skeleton.PGC ABC3.Found.GaloisRep
    open scoped NNReal Valued

の 3 行を付ける（`𝒪[K.carrier]` / `𝓀[K.carrier]` は `open scoped Valued` の記法）。
★対象ファイルの `open` 行をまず `grep -n '^open' <file>` で写すのが確実。

## #82 `PowerSeries` の補題が `MvPowerSeries.map` / `MvPowerSeries.expand` で書かれていて `rw` が刺さらない（2026-09-06、Λ6 DworkThetaStep2）

**症状**: `PowerSeries.map_subst` / `PowerSeries.expand_subst` /
`MvPowerSeries.map_iterateFrobenius_expand` を `rw` すると

    Tactic `rewrite` failed: Did not find an occurrence of the pattern
      PowerSeries.subst ((PowerSeries.expand (pp ^ n) hqne) Z) ?m
    in the target expression
      PowerSeries.subst ((MvPowerSeries.expand (pp ^ n) hqne) Z) P

が出る。目にはまったく同じ式に見えるのに一致しない。

**原因**: `PowerSeries R := MvPowerSeries Unit R` で、`PowerSeries.map` は
`MvPowerSeries.map` の、`PowerSeries.expand` は `MvPowerSeries.expand` の
**別名（定義そのもの）**。定義的には等しいが `rw` の keyed matching は
ヘッド定数が違うと当たらない。mathlib の `PowerSeries.*` 補題は
statement の一部を `MvPowerSeries.*` のまま書いているものがある。

**直し方**: 使いたい形を `have` で 1 度だけ言い直す（型を明示すれば defeq で通る）。

    have h1 : PowerSeries.map (iterateFrobenius R pp n)
        (PowerSeries.expand (pp ^ n) hqne W) = W ^ (pp ^ n) :=
      MvPowerSeries.map_iterateFrobenius_expand (σ := Unit) (R := R) pp hp W n

`map_subst` は 1 変数版を 1 本作って以後それだけ使うのが安い。

    theorem map_subst_powerSeries (φ : B →+* C) {a f : PowerSeries B}
        (ha : PowerSeries.HasSubst a) :
        PowerSeries.map φ (PowerSeries.subst a f)
          = PowerSeries.subst (PowerSeries.map φ a) (PowerSeries.map φ f) :=
      PowerSeries.map_subst ha f

## #83 `set x := (PowerSeries の式)` のあと `x.foo` が `Function.foo` を探しにいく（2026-09-06、Λ6 DworkThetaStep2）

**症状**:

    set θ0 := PowerSeries.subst θ M with hθ0def
    refine ⟨θ0, θ0.substInvOfIsUnit hθ01u, …⟩

で

    Invalid field `substInvOfIsUnit`: The environment does not contain
    `Function.substInvOfIsUnit` … from an expression θ0 of type
    (Unit →₀ ℕ) → ↥(unramifiedCompletionInt K)

**原因**: `PowerSeries.subst` の戻り値の型は `MvPowerSeries τ S` と書かれており、
`set` は**推論された型そのまま**（`MvPowerSeries Unit B`、さらに展開すると関数型）
で局所定義を作る。dot notation は型のヘッドを見るので `PowerSeries` に届かない。

**直し方**: どちらか。

* `set θ0 : PowerSeries B := … with h` と**型を明示**する（これで `set` は
  `PowerSeries` のまま持つ）。
* dot notation をやめて `PowerSeries.substInvOfIsUnit θ0 h` と完全修飾で書く。

★型注釈を付けても、`PowerSeries.constantCoeff_subst_eq_zero` のように
結論が `MvPowerSeries.constantCoeff` で書かれた補題は #82 の問題が残る。

## #84 `τ * σ * τ⁻¹ • y` は `τ * σ * (τ⁻¹ • y)` に読まれる（2026-09-06、Y1 LowerRamificationGroup）

**症状**: 正規部分群の共役条件を書こうとして

    have h1 : τ * σ * τ⁻¹ • y - y = τ • (σ • (τ⁻¹ • y) - τ⁻¹ • y) := …

で `failed to synthesize HMul G B`。`G` の元と環の元を掛けようとしている。

**原因**: `•` は `*` より**結合が強い**。`τ * σ * τ⁻¹ • y` は
`τ * σ * (τ⁻¹ • y)` と解析される。

**直し方**: `(τ * σ * τ⁻¹) • y` と括る。`Subgroup.Normal` の `conj_mem` を
`Ideal.inertia` で書くときに毎回踏む。

## #85 statement に `IsDiscreteValuationRing.addVal B` を出すと `haveI` では遅い（2026-09-06、Y1）

**症状**:

    theorem foo … : … (n : ℕ∞) < IsDiscreteValuationRing.addVal (adjoinIntegers K x) z := by
      haveI := isDiscreteValuationRing_adjoinIntegers K x   -- ← 遅い
      …

で `failed to synthesize IsDiscreteValuationRing ↥(adjoinIntegers K x)`。

**原因**: インスタンスが要るのは**証明中ではなく statement の elaboration 時**。
`haveI` は本体に入ってからしか効かない。

**直し方**: ファイル（節）の先頭で

    attribute [local instance] isDiscreteValuationRing_adjoinIntegers

と入れる。ABC3 では DVR 性が `instance` ではなく `theorem` で置かれている
（`isDiscreteValuationRing_adjoinIntegers` / `isDiscreteValuationRing_carrierIntegers` /
`module_finite_adjoinIntegers`）ので、statement に出す節でだけ local instance にする。

## #86 `Submodule.smul_induction_on` の `a • b` は `smul_eq_mul` を先に要る（2026-09-06、Y1）

**症状**: `x ∈ I * J`（`Ideal`）に `Submodule.smul_induction_on` を当てると
ゴールが `σ • a • b ∈ …` になり、`rw [smul_mul']` が
「`?a • (?b₁ * ?b₂)` が見つからない」で落ちる。

**原因**: `Ideal` の積は `Submodule` の `smul` で定義されているので、
分解して出てくるのは `a • b`（環自身へのスカラー倍）であって `a * b` ではない。

**直し方**: `rw [smul_eq_mul, smul_mul']` と、**先に `smul_eq_mul` で `*` に直す**。

## #87 `Ideal.span_singleton_pow` の向き（2026-09-06、Y1）

`Ideal.span_singleton_pow : Ideal.span {a} ^ n = Ideal.span {a ^ n}`。
`𝔪 = span {α}` を `rw` したあと `𝔪 ^ k` は `span {α} ^ k` になるので、
`Ideal.mem_span_singleton'` に持ち込むには `rw [Ideal.span_singleton_pow]`（← を付けない）。

## #88 `congrArg Multiplicative.ofAdd ?_` は defeq をすり抜けて元のゴールを返す（2026-09-06、Y3）

**症状**: `M →* Multiplicative R` の `map_mul'` や `MonoidHom.ext` の後で
`refine congrArg Multiplicative.ofAdd ?_` と書くと**成功したように見えて**、
残ったゴールが元のまま（あるいは `toAdd (ofAdd a + ofAdd b)` という珍妙な形）になり、
次の `rw` が「パターンが無い」で落ちる。`Add (Multiplicative R)` が無い、という
無関係なエラーが出ることもある。

**原因**: `Multiplicative R` は `R` と defeq なので、`congrArg ⇑Multiplicative.ofAdd`
の `?a`/`?b` がゴールの左右にそのまま当たってしまい、`ofAdd` を剥がせない。

**直し方**: 剥がした後の等式を `have h : … = … := by …` で**先に作ってから**
`exact congrArg Multiplicative.ofAdd h`。`1` を消したいだけなら `ofAdd_eq_one`
（`ofAdd x = 1 ↔ x = 0`）を `rw`/`simp only` する。

## #89 `isUnit_of_mul_eq_one` は `IsUnit.of_mul_eq_one`、`add_right_eq_self` は `add_eq_left`（2026-09-06、Y3）

現行 mathlib での名前。`IsUnit.of_mul_eq_one (b) (h : a * b = 1) : IsUnit a`
（`[IsDedekindFiniteMonoid M]` が要るが可換なら自動）。

## #90 `MulSemiringAction` の `σ • x ^ n` は `smul_pow'`（2026-09-06、Y3）

`smul_pow' : r • x ^ n = (r • x) ^ n`（`MulDistribMulAction`、`@[simp]`）。
`smul_pow : (r • x) ^ n = r ^ n • x ^ n` は**別物**（`Monoid` へのスカラー倍）で、
`rw [smul_pow]` は `(?r • ?x) ^ ?n` を探して落ちる。

## #91 `rw [← mul_smul]` は 1 回だと片側しか畳まず `congr 1` が誤分割する（2026-09-06、Y4）

`(τ * σ * τ⁻¹) • (τ • α) = τ • (σ • α)` を出したいとき

```lean
rw [← mul_smul]; congr 1; group   -- ✗
```

は左辺だけ `(τ*σ*τ⁻¹*τ) • α` に畳み、右辺は `τ • σ • α` のままなので
`congr 1` が `τ * σ = τ` と `α = σ • α` という**偽のゴール 2 本**に割る。

**直し方**: 両側を畳んでから割る。

```lean
rw [← mul_smul, ← mul_smul]; congr 1; group   -- ○
```

## #92 群 `G` が結論に現れない補題は `(G := G)` を明示しないと instance が stuck（2026-09-06、Y4）

`exists_lowerRamificationGroup_eq_bot (A := A) hadj` のように、結論
`∃ N, ∀ n, N ≤ n → lowerRamificationGroup B G n = ⊥` を
`obtain` の右辺でしか使わない形だと、`G` が未決定のまま
`[FaithfulSMul ?m B]` の探索に入り

```
typeclass instance problem is stuck
  FaithfulSMul ?m.49 B
```

で止まる。`hadj : Algebra.adjoin A {α} = ⊤` は `G` を決めない。
**直し方**: `(A := A) (G := G)` と**両方**明示する。

## #93 `[IsDomain B] [IsLocalRing B] [IsDiscreteValuationRing B]` を並べると `linter.overlappingInstances` が鳴る（2026-09-06、Y5）

DVR の補題を書くときに「必要なものを全部並べる」と

```
⚠️ `[IsDomain B]`, `[IsLocalRing B]`, and `[IsDiscreteValuationRing B]` each imply `[Nontrivial B]`.
💡️ Of these, `[IsLocalRing B]` may be removed.
```

が出る。`IsDiscreteValuationRing` は `IsLocalRing` を含むので **`[IsLocalRing B]` を落とす**。
★落としても `maximalIdeal B` は書ける（DVR インスタンスから出る）。
★具体層で `attribute [local instance] isDiscreteValuationRing_adjoinIntegers` を入れて
`lowerRamificationGroupAdjoin`（`instIsLocalRingAdjoinIntegers` で作った定義）に
`exact` する経路も通る——`LowerRamificationGroup.lean` の
`mem_lowerRamificationGroupAdjoin_iff_lt_addVal` が先例。

## #94 `rw [mul_comm]` は `∣` の**左辺**に当たる（2026-09-06、Y5）

ゴール `α * π ∣ ↑m * α` で右辺だけ入れ替えたくて `rw [mul_comm]` と撃つと
`π * α ∣ ↑m * α` になり、`mul_dvd_mul_left α hpi : α * π ∣ α * ↑m` が当たらない。
**直し方**: `rw [mul_comm (m : B) α]` と**引数を書く**。

## #95 `Function.update_noteq` / `update_same` は `update_of_ne` / `update_self`、`σ • ∏` は `Finset.smul_prod'`（2026-09-07、Y6）

桁展開の帰納で `Function.update a M d` を作るとき、旧名 `Function.update_noteq` /
`Function.update_same` は**もう無い**（`Unknown constant`）。現行名は
`Function.update_of_ne : a ≠ b → Function.update f b v a = f a` と
`Function.update_self`。
同じ回で `σ • ∏ i ∈ s, f i = ∏ i ∈ s, σ • f i` を `Finset.prod_map_smul` で探して
外した。正しくは **`Finset.smul_prod'`**（`MulDistribMulAction` 版。`'` が付く方）。
**直し方**: 名前が出てこない補題は `.cache/mathlib-index.txt` を結論の形で引く。
「作用が積を保つ」は `smul_prod'`、「作用が和を保つ」は `map_sum`
（`MulSemiringAction.toRingHom G B σ` に書き換えてから）。

## #96 `set` で束縛した局所定義は `rw` のパターンに当たらない —— `have h2 : <展開形> = _ := h` を挟む（2026-09-07、Λ6a′）

`set w' : adjoinIntegers K w := ⟨⟨w, hwmem⟩, _⟩` としたあと、補題から得た
`hE : ... = ↑↑w'` でゴール `w ∈ K⟮u⟯` を `rw [← hE]` しようとすると
`Did not find an occurrence of the pattern ↑↑w'` で落ちる（ゴールには
`↑↑w'` ではなく既に簡約された `w` が現れているため）。
**直し方**: `have h2 : (↑↑(...) : K.closure) = w := hE` と**欲しい形で型を書き直して**
（`↑↑w'` と `w` は defeq なので `exact hE` で通る）から `rw [← h2]`。
★同じ回で `congr 1` も外した——2 段の coe（`adjoinIntegers → K⟮x⟯ → K.closure`）に
撃つと 1 段しか剥がれず向きも反転する。`congr` ではなく
「等式を先にゴールへ `rw` して `exact h.symm`」が速い。

## #97 `K.carrier` から `K` は逆算できない —— `(K := …)` を名前付き引数で先に固定する（2026-09-07、Λ7 D24 第 1 段）

`h : ∀ {K K' : PAdicLocalField p}, (K.carrier ≃ₐ[ℚ_[p]] K'.carrier) → …` に
`twistedAlgEquiv p : TwistedQp p ≃ₐ[ℚ_[p]] ℚ_[p]` をそのまま当てると
`Application type mismatch … expected PAdicLocalField.carrier ?m ≃ₐ[ℚ_[p]] PAdicLocalField.carrier ?m'`
で落ちる（構造体の射影 `carrier` を逆向きに解くのは高階単一化で、`?m` が決まらない）。
**直し方**: `h (K := twistedField p) (K' := selfField p) (twistedAlgEquiv p)` と
名前付き引数で `K`・`K'` を先に固定する。固定してしまえば本体は defeq で通る
（`(selfField p).carrier` は `exact` なら `ℚ_[p]` と合う、の項と同じ事情）。

## #98 MCP REPL の基準環境は**共有**である —— 並行 agent に import を差し替えられる（2026-09-07、Y7a）

`lean_start(["ABC3.Found.PGC.UniformizerExpansion"])` が **10.2 秒で「成功」**と返ったのに、
その後の `lean_check` で `ABC3.Found.PGC.ramIndex` が `Unknown identifier` になった。
`lean_status` を見ると imports が **`ABC3.Check.PGC.Prop12Degenerate, …`**（別 agent のもの）で、
基準環境が横から差し替わっていた。**起動が速すぎる（90 秒でなく 10 秒）ときは疑うこと。**
**直し方**: `lean_start` を撃ち直すと相手の環境を壊して往復戦争になる（役割定義の 313 回問題）。
**`node tools/leanfile.mjs <path>` に切り替える**（olean を書かないので並行安全、12 秒/往復）。
★実測: Y7a は §1 の 4 補題だけ MCP で通し、以降は leanfile.mjs で 3 往復。ファイル 455 行は初回で通った。

## #99 `rw` の末尾 `rfl` は reducible 止まりなので、部分体の coe のノルム一致を閉じられない（2026-09-07、M3）

`E : IntermediateField K L` の元 `x` について `‖x‖ = ‖(x : L)‖` は **`rfl` で通る**
（`↥E` のノルムは `L` のノルムの制限そのもの）。ところが

    rw [foo, bar]          -- ゴールが `‖↑x‖ = ‖x‖` になって止まる

`rw` が最後に試す `rfl` は **reducible 透明度**なので、`SubfieldClass` 経由のノルムの
定義を開けず閉じられない。エラーは `unsolved goals ⊢ ‖↑x‖ = ‖x‖` という、
一見「同じ式なのに」に見える形で出る。
★直し方は**名前付き補題を 1 本置いて `rw` の鎖に入れる**（`rfl` を単独のタクティクとして
書き足してもよいが、同じ形が何度も出るので補題にする方が安い）:

    theorem norm_coe_sub (x : ↥E) : ‖(x : L)‖ = ‖x‖ := rfl   -- 0.03 秒
    ... rw [foo, bar, norm_coe_sub]

## #100 `lean_start` が 10 秒で返っても、`lean_status` の imports が自分のものなら正常（2026-09-07、M3）

#98 は「10 秒で返ったら他 agent の環境を掴んでいる」と書いたが、**それだけでは判定にならない**。
実測（M3、並行 2 体）: `lean_start(["…UnramifiedCompletion","…AdjoinIntegers"])` が **10.6 秒**で
返り、`imports:` は**要求どおり**、`#check` も自分の宣言が全部見えた（olean が OS の
ページキャッシュに乗っていた）。★**判定は所要時間ではなく `lean_status` の imports 一致で行う。**
一致していれば撃ち直さない（撃ち直しがそのまま費用になる）。


## #101 `rw [iff補題]` はゴールが `≠` だと当たらない —— `intro` を先に打つ（2026-09-07、Y7b）

`ramIndex_pow_pow_eq_top_iff : ramIndex π (σ ^ p ^ k) = ⊤ ↔ m ≤ k` を、ゴール
`ramIndex π (σ ^ p ^ k) ≠ ⊤` に `rw` で当てると

    error: Tactic `rewrite` failed: Did not find an occurrence of the pattern

になる。`Ne a b` は `a = b → False` で、`rw` は `Ne` を展開しないまま
`a = b` を探すため（`ne_eq` を `simp only` で剥がすか、下のように `intro` する）。

    -- ✗
    refine ENat.coe_toNat ?_
    rw [ramIndex_pow_pow_eq_top_iff (A := A) hp.one_lt hadj hord k]
    omega
    -- ✓
    refine ENat.coe_toNat ?_
    intro hc
    have h2 := (ramIndex_pow_pow_eq_top_iff (A := A) hp.one_lt hadj hord k).1 hc
    omega

★同じ形は `x ∉ S`（= `x ∈ S → False`）でも起きる。**`≠` / `∉` のゴールは
`rw` ではなく `intro` して `.1` / `.2` を当てる**のが速い。

## #102 `ℕ∞` の引き算は切り詰めるので「合同」を書くと空虚に真になる（2026-09-07、Y7b）

`a < b` なら `ℕ∞` で `a - b = 0`（`tsub_eq_zero_of_le`）。したがって
`(p ^ j : ℕ∞) ∣ (i_j - i_{j+1})` は **`simp` で通る内容ゼロの定理**である
（実測: `example (a b c : ℕ∞) (h : a < b) : c ∣ a - b` が 2 行で閉じる）。
★`ℕ∞` 値の量の「合同」を主張するときは、**有限代表 `c d : ℕ` を仮定
（`f j = (c : ℕ∞)`）で取り出し、`ℤ` の中の `(d : ℤ) - (c : ℤ)` で書く**こと。
`⊤` の側は「`f (j+1) = ⊤ ∨ ∃ c d, …`」の場合分けで明示する。

## #103 `omega` は積を**向きごと**に別の原子として数える（2026-09-07、Y8）

`obtain ⟨c, hc⟩ := (h : p ^ (k+1) ∣ f (k+1) - f k)` が返す `hc` は
`f (k+1) - f k = p ^ (k+1) * c`（**割る数が左**）。一方ゴールは
`f (k+1) = f k + c * p ^ (k+1)`（**係数が左**）と書きたい。この 2 つを一緒に `omega` へ渡すと

    omega could not prove the goal:
     g := ↑(p ^ (k + 1)) * ↑c
     h := ↑c * ↑(p ^ (k + 1))

と、**同じ積を 2 つの原子 `g` `h` として並べたまま**落ちる（`g - h ≤ -1` が反例に出る）。
`omega` は線形算術なので、非線形な積は原子として**構文で**同一視する。

    -- ✓ 先に向きを揃えてから omega に渡す
    have hc' : f (k + 1) - f k = c * p ^ (k + 1) := by rw [hc, mul_comm]
    omega

★`ring_nf at hc ⊢` でもよいが、`have` で**欲しい向きの等式を 1 本置く**方が読みやすい。

## #104 仮定から述語を逆算させない —— `(P := fun k => …)` を明示する（2026-09-07、Y8）

`theorem exists_index_boundary {P : ℕ → Prop} (h0 : P 0) : ∀ M, ¬ P M → …` に
`exists_index_boundary h0 m hm` と当てると

    Application type mismatch: hm has type ¬ramIndex π (σ ^ p ^ m) ≤ ↑n
    but is expected to have type ¬?m.312 m

になる。`h0 : ramIndex π (σ ^ p ^ 0) ≤ ↑n` から `P` を作るのは高階の逆算で、
Lean は `P := fun _ => (その命題)` の方を先に選んでしまう。

    -- ✓
    exists_index_boundary (P := fun k => ramIndex π (σ ^ p ^ k) ≤ (n : ℕ∞)) h0 m hm

★#92（群 `G` が結論に出ない補題）と同じ「明示に固定する」対処。**述語を引数に取る補題は
呼ぶ側が必ず `(P := …)` を書く**と決めておくと往復が 1 回減る。
★ついでに `Subgroup.eq_top_iff'` は部分群を**明示引数**で取るので `Subgroup.eq_top_iff'.2 hg`
は `Invalid projection: Projections cannot be used on functions` になる。`(Subgroup.eq_top_iff' _).2 hg` と書く。

## #105 statement に instance が要るなら `haveI` では間に合わない —— `attribute [local instance]`（2026-09-07、Λ#7）

`PowerSeries.aeval hz θ` を **statement に**書くと、その型に `IsLinearTopology S S` や
`ContinuousSMul A S` が要る。`ClosureCompletion.lean` はこれらを `theorem` として持っていて
（`instance` ではない）、`def` の中では `haveI := isLinearTopology_closureCompletionInt K` と
借りている。同じ手を **theorem の statement** でやると

    theorem foo … : ‖PowerSeries.aeval hz θ‖ ≤ ‖z‖ := by
      haveI := isLinearTopology_closureCompletionInt K   -- ← 遅い
    → failed to synthesize instance  IsLinearTopology ↥(closureCompletionInt K) …

になる（statement の elaborate は `haveI` より前）。#85「`haveI` では遅い」の instance 版。

    -- ✓ ファイル冒頭（namespace の中）で 1 回
    attribute [local instance] isLinearTopology_closureCompletionInt continuousSMul_closureCompletionInt

`local` なのでファイルの外へは漏れない（他ファイルの instance 探索を汚さない）。
★`attribute` の直前に `/-- … -/` を置くと `unexpected token 'attribute'; expected 'lemma'`。
説明は `/-! … -/` で書く。

## #106 `Polynomial.coe_map` は無い —— `Polynomial.polynomial_map_coe`（2026-09-07、Λ#7）

`PowerSeries.map φ ↑P = ↑(P.map φ)` が欲しいときの名前は
`Polynomial.polynomial_map_coe : Polynomial.map φ f = PowerSeries.map φ ↑f`
（`RingTheory/PowerSeries/Basic.lean`）。`← Polynomial.polynomial_map_coe` で使う。
`Polynomial.coe_map` は `Unknown constant`。
★ついでに `Multiset.prod_eq_zero_iff.mp h` が返すのは **`0 ∈ s` という所属**であって
`⟨w, hw, hw0⟩` に分解できるタプルではない（`rcases failed: … is not an inductive datatype`）。
`Multiset.mem_map.mp (Multiset.prod_eq_zero_iff.mp h)` と繋ぐ。

## #105 `Algebra.adjoin_singleton_eq_range_aeval` の `obtain` は `.toRingHom` を被せて返す（2026-09-07、Y9）

`rw [Algebra.adjoin_singleton_eq_range_aeval] at hy; obtain ⟨q, hq⟩ := hy` の `hq` は

    hq : (aeval t).toRingHom q = y

であって `(aeval t) q = y` **ではない**（`AlgHom.range` の membership が
`RingHom.range` 経由で展開される）。そのため直後の

    rw [← hq, aeval_def, eval₂_eq_eval_map]   -- ✗

は `Did not find an occurrence of the pattern (aeval ?x) ?p` で落ちる。

    -- ✓ defeq で型を付け替える 1 行を挟む
    have hq' : (aeval t) q = y := hq
    rw [← hq', aeval_def, eval₂_eq_eval_map]

★同じ形は `AlgHom.range` / `RingHom.range` / `MonoidHom.mrange` の `obtain` 全般で起きる。
**`rw` が「無いはずのない項が見つからない」と言ったら、まず `have` で defeq に付け替える。**
`simp only [AlgHom.toRingHom_eq_coe, RingHom.coe_coe]` でも直るが、`have` の方が 1 行短く速い。

## #107 `open ... in` は docstring の**前**に置く（2026-09-07、Y9）

```lean
-- x  unexpected token 'open'; expected 'lemma'
/-- 説明 -/
open scoped PowerSeries.WithPiTopology in
theorem foo ...

-- o
open scoped PowerSeries.WithPiTopology in
/-- 説明 -/
theorem foo ...
```

★#105 の `attribute` と同じ形の罠（`/-- … -/` の直後に来られるのは宣言だけ）。
`lean_check` では docstring を省いて通していたので**ファイルに書いた瞬間に初めて出る**
——docstring を足したら `lake build` の前にもう一度 `lean_check` に通すこと。

## #108 完備体上の「有限次元なら閉」は `NormedSpace` を要求しない（2026-09-07、Y9）

`Submodule.closed_of_finiteDimensional` の仮定は

    [NontriviallyNormedField 𝕜] [CompleteSpace 𝕜] [AddCommGroup E] [TopologicalSpace E]
    [IsTopologicalAddGroup E] [Module 𝕜 E] [ContinuousSMul 𝕜 E] [T2Space E]

で、**`NormedSpace 𝕜 E` も `E` のノルムも要らない**。完備化 `ℂ_K` の上で
「`K̂^{ur}` 上有限次元の中間体は閉」を出すのに `NormedSpace K̂^{ur} ℂ_K` を
組もうとすると `‖a • x‖ = ‖a‖‖x‖` の証明で時間を溶かす——**組まなくてよい**。
必要なのは `Algebra`（RingHom から `.toAlgebra`）と `ContinuousSMul`（`continuous_mul`
と埋め込みの連続性の合成）の 2 つだけ。中間体版は 8 行:

    theorem isClosed_intermediateField_of_finiteDimensional (E : IntermediateField 𝕜 L)
        [FiniteDimensional 𝕜 E] : IsClosed (E : Set L) := by
      haveI : FiniteDimensional 𝕜 ↥(Subalgebra.toSubmodule E.toSubalgebra) :=
        inferInstanceAs (FiniteDimensional 𝕜 E)
      have h := Submodule.closed_of_finiteDimensional (Subalgebra.toSubmodule E.toSubalgebra)
      rwa [Subalgebra.coe_toSubmodule] at h

★`Algebra.fg_adjoin_of_finite` は **`Algebra.` が付かない**（根名前空間の
`fg_adjoin_of_finite`）。中間体なら `IntermediateField.finiteDimensional_adjoin`
（`[Finite S]` が instance 引数なので `haveI := hS.to_subtype` を先に置く）。

## #107 2026-09-07 に名前が動いていた 5 つ（Y10）

`lean_check` を使わず `node tools/leanfile.mjs` だけで 8 往復した際に踏んだもの。
**どれも「昔の名前を書くと `Unknown identifier` / deprecated 警告」**で、直し方は名前の置換だけ。

| 書きたかったもの | 通る名前 |
| --- | --- |
| `hirr.not_unit`（`Irreducible` の構造体フィールド） | **`hirr.not_isUnit`** |
| `rcases le_or_lt a b with h \| h` | `le_or_lt` は無い。**`by_cases h : a ≤ b` + `Nat.not_le.mp h`** |
| `isUnit_of_mul_eq_one a b h` | **`IsUnit.of_mul_eq_one b h`**（第 1 引数が消えた） |
| `Algebra.algebraMap_mem S r` | **`Subalgebra.algebraMap_mem S r`** |
| `Polynomial.eval_finset_sum` / `ENat.one_le_iff_ne_zero` / `mul_le_mul_left'` / `push_neg` | `Polynomial.eval_finsetSum` / `Order.one_le_iff_ne_zero` / `mul_le_mul_right` / `push Not` |

★**`IsIntegral R x` を無名構成子で開くと、残る goal は `aeval` ではなく `eval₂` である。**

    refine ⟨q, hqm, ?_⟩
    -- ⊢ eval₂ (algebraMap C B) π q = 0     ← `(aeval π) q = 0` ではない
    rw [Polynomial.aeval_def]            -- ✗ Did not find an occurrence of the pattern (aeval π) q
    rw [Polynomial.eval₂_eq_eval_map]    -- ✓ そのまま `map` に移れる

同じ「`rw` が無いはずのない項を見つけられない」形は #105 と同種。**まず goal を見る。**

## #109 ℤ 上の両方向の帰納法まわりで踏んだ 4 つ（2026-09-07、Y4b / Lemma 4.5）

`Yoshida08 Lemma 4.5`（`θ^{(j)}/θ = π′_j/π_j`、`j : ℤ`）を書くときに踏んだもの。
**どれも「ℤ を負の側まで通す」ときにだけ出る。**

| 書きたかったもの | 通る形 |
| --- | --- |
| `rcases lt_or_le j 0 with hj \| hj` | `lt_or_le` は無い。**`rcases (by omega : j < 0 ∨ 0 ≤ j) with hj \| hj`**（omega が選言そのものを証明する。ℤ/ℕ ならこれが一番壊れない） |
| `induction j using Int.induction_on with \| hz \| hp \| hn` | 場合の名前は **`zero` / `succ` / `pred`**（`hz`/`hp`/`hn` は古い。`Invalid alternative name` と出る） |
| `Units.val_prod` | 無い。**`map_prod (Units.coeHom L) f s`** |
| `(-(n:ℤ)).toNat = 0` / `(-(-(n:ℤ))).toNat = n` | **`by omega`** が両方通す（`Int.toNat_*` の名前を探さない） |

★**`map_prod (Units.coeHom L) …` を `rw` に渡すと必ず失敗する。**
goal には `↑(∏ …)`（`Units.val` の coercion）が出ていて `(Units.coeHom L) (∏ …)` は
**構文的に一致しない**（defeq ではある）。`refine (map_prod …).trans ?_` にすれば通る:

    rw [uniformizerZ, zpowProd_natCast, natProd]
    refine (map_prod (Units.coeHom L) (fun t => ((unitsRingAutHom L ϕ) ^ t) π)
      (Finset.range n)).trans ?_
    exact Finset.prod_congr rfl fun t _ => coe_unitsRingAutHom_pow ϕ t π

★**`ring` / `linear_combination` は defeq を見ない。**
`↑(σ u)`（`Units.val` を通した像）と `ϕ ↑u` は defeq だが `ring` には別の原子に見え、
`ring failed, ring expressions not equal` になる。**`rfl` の `have` を 1 本挟んで
`rw` で構文を揃えてから** `linear_combination` を呼ぶ:

    have hcoe : ((((unitsRingAutHom L ϕ) x) : Lˣ) : L) = ϕ ((x : Lˣ) : L) := rfl
    rw [hcoe] at hsub
    linear_combination -hsub

★**環自己同型を `ℤ` 乗したいときは `MonoidHom` を 1 本作る。**
`RingAut R →* MulAut Rˣ`（`toFun ϕ := Units.mapEquiv ϕ.toMulEquiv`、`map_one'`/`map_mul'`
はどちらも `by ext u; rfl`）を作っておけば `map_zpow` が使えて
`((σ^j) u : R) = (ϕ^j) (u : R)` が `rw [← map_zpow …]; rfl` の 1 行で出る。
**符号の場合分けが要らなくなる**のが効き目。

## #110 `Subgroup G` を `Set G` に落とすとき `G` が推論できない（2026-09-07、Y11 / Prop 6.9）

`Subgroup` の carrier を `Set` として書く statement

    (ramificationGroupReal π' n : Set G) * (H : Set G) = …

は **`Type mismatch: … has type Subgroup ?m.130 but is expected to have type Set G`**
で落ちる。`ramificationGroupReal α n : Subgroup G` の `G` は
`[MulSemiringAction G B]` からしか決まらないので、`: Set G` という型上昇の**外側**の
`G` が中の `?m` に伝わらない。★`SetLike` の coercion は「先に `Subgroup ?G` を
elaborate してから `Set ?G` に落とす」ので、`Set G` の `G` は unifier に見えていない。

通る形（2 段に分けて `G` を明示する）:

    ((ramificationGroupReal (G := G) π' n : Subgroup G) : Set G) * (H : Set G)

★`(H : Set G)` は `H : Subgroup G` が既に explicit な変数なので落ちない。
落ちるのは**返り値の型引数が implicit な関数**の像だけである。

## #111 `field_simp [f]` は `def f` を展開しない（2026-09-07、Y11）

`herbrandPhi α H n := -1 + herbrandSum α H n / (Nat.card H : ℝ)` に対し

    field_simp [herbrandPhi]

は `herbrandPhi` を展開せずに `unsolved goals` を残す（simp 引数に def 名を渡しても
equation lemma が使われない場面がある）。**`simp only [herbrandPhi]` を先に打ってから
`field_simp`**、さらに `field_simp` が `-c + S + c = S` の形で止まるので **最後に `ring`**。

    simp only [herbrandPhi]
    field_simp
    ring

★同じ落とし穴は `rw [herbrandSum]` にもある（`unfold herbrandSum` なら通る）。

## #112 `push_neg` が 2026-09-07 に deprecated になった（Y11）

`by_contra hc; push_neg at hc` は **warning: `push_neg` has been deprecated.
Prefer using `push Not` instead.** を出す。`ℕ` の `¬ (a < b)` を返すだけなら
**`Nat.not_lt.1 hc`**（一般には `not_lt.1`）で置き換えるのが一番安い。

★`lake build` は warning では落ちないが、`leanfile.mjs` の出力が warning で
埋まって本物のエラーが見えにくくなる。

## #113 ℤ で添字づけた `Σ` 型（非交和）を作るときの 4 つの穴（Y-Prop4.7）

**(a) `Equiv.sigmaFiberEquiv` と `Equiv.sigmaCongrLeft'` の合成は `rfl` にならない。**
`(Equiv.sigmaFiberEquiv w).symm.trans (Equiv.sigmaCongrLeft' (Equiv.neg ℤ))` で
`X ≃ Σ j, {x // w x = -j}` は作れるが、`sigmaCongrLeft'` が `Eq.mpr` を挟むので
`apply` 補題が `rfl` で証明できない（`Not a definitional equality`）。
**`Equiv` を手で組む方が安い**（`toFun`/`invFun` を明示すれば `_apply` が全部 `rfl`）。
`right_inv` の依存等式は **`obtain ⟨x, hx⟩ := y; have hj : j = -w x := by omega; subst hj; rfl`**
——添字 `j` を `subst` で潰せば証明無関係で `rfl` になる。

**(b) `Function.Bijective.1` は β 簡約した等式を受け取らない。**
`hbij : Function.Bijective (fun g : S => (⟨ρ g, _⟩ : T))` に
`Subtype.ext h : (⟨ρ a, _⟩ : T) = ⟨ρ b, _⟩` を渡すと
`Application type mismatch`（期待は `(fun g => …) ?a = (fun g => …) ?b`）。
**`have h' : (fun g : S => …) ⟨a, ha⟩ = (fun g : S => …) ⟨b, hb⟩ := Subtype.ext hab`**
と、ラムダを書いたままの型で `have` を立ててから渡す。

**(c) `MonoidHom.ker` の標的は乗法群でなければならない。**
`v : Q →* ℤ` の `v.ker` は通らない（ℤ は乗法モノイドで群でない）。
付値は **`Q →* Multiplicative ℤ`** で持つこと。
`{x // Multiplicative.toAdd (w x) = -j}` と `{y // w y = Multiplicative.ofAdd (-j)}` は
**defeq なので `Equiv.trans` がそのまま通る**（型合わせの `cast` は要らない）。

**(d) `MulEquiv.ofBijective` は noncomputable。**
`def foo … : G ≃* H := MulEquiv.ofBijective ρ h` は
`failed to compile definition, consider marking it as 'noncomputable'` で落ちる。
★**全単射性の `theorem` と `≃*` の `noncomputable def` を分けて書く**と、
下流が `Bijective` だけ欲しいときに `noncomputable` が伝染しない。

## #114 束縛子の型を書かないと `(τ : G)` が「coe」ではなく「型注釈」に読まれる（Y12・Lemma 6.10）

`H : Subgroup G` の上で立てた補題に渡す関数を

```lean
rw [phiOf_natCast_of_pos _ (fun τ => pos_ramIndex huni (τ : G)) n]
```

と書くと、`τ` の型がまだメタ変数なので **`(τ : G)` が型注釈として働き `τ : G` に固定される**。
その結果 `phiOf_natCast_of_pos` の `{S}` が `S := G` に解かれ、
**`failed to synthesize Fintype G`**（本来は `Fintype ↥H` で足りるはず）が出る。
エラーが「`Fintype G` が無い」と言うので**インスタンスの問題に見えるが、原因は elaboration の順序**である。

★直し方: **束縛子に型を書く**。

```lean
rw [phiOf_natCast_of_pos _ (fun τ : H => pos_ramIndex huni (τ : G)) n]
```

こうすると `(τ : G)` は `↥H → G` の coe として読まれ、`S := ↥H` に解ける。
★同じ形は `Finset.sum_congr rfl (fun τ => …)` / `Fintype.sum_equiv … (fun τ => …)` でも起きる。
**部分群の上の和を扱うラムダは、常に `fun τ : H =>` と書くこと。**

## #115 `ℕ` の区間和で踏んだ 3 つ（Y13・Corollary 6.3 / Hasse-Arf 受け渡し）

**(a) `Nat.Ico_succ_right` は無い。** `Finset.Icc 1 n` と `Finset.Ico 1 (n+1)` を行き来したいとき
`rw [← Nat.Ico_succ_right]` は `Unknown constant` で落ちる（2026-09-07 時点）。
在庫は `Finset.Ico_succ_right_eq_Icc` / `Order.Ico_succ_right` だが、どちらも `Order.succ` 版で
`ℕ` に当てるのに一手要る。★**その場で作る方が速い**:

```lean
have h : Finset.Icc 1 (e * m) = Finset.Ico 1 (e * m + 1) := by
  ext x; simp only [Finset.mem_Icc, Finset.mem_Ico]; omega
```

**(b) `Nat.mul_le_mul_left e (Nat.le_succ m)` は `e * m.succ` を作り、`omega` が
`e * (m + 1)` と別物として扱う。** 帰納法の `succ` 分岐で

```lean
have hle : e * m + 1 ≤ e * (m + 1) + 1 := by
  have := Nat.mul_le_mul_left e (Nat.le_succ m); omega   -- ✗ 反例を出して落ちる
```

とすると `a := ↑e * ↑m.succ` と `c := ↑e * ↑(m + 1)` が**別の原子**になり `omega` が
「反例があるかもしれない」と言う（掛け算の中は `omega` が正規化しない）。★直し方は
**`Nat.mul_succ` で先に展開する**:

```lean
have hle : e * m + 1 ≤ e * (m + 1) + 1 := by rw [Nat.mul_succ]; omega   -- ○
```

同じ理由で `e * (m+1) ≤ e * k` から `e*m + e ≤ e*k` を出すのも `rw [Nat.mul_succ] at h` が要る。

**(c) `exists_prime_orderOf_dvd_card`（Cauchy）は `Fintype` を要求し、`Nat.card` を受けない。**
`[Finite Q]` と `p ∣ Nat.card Q` から使うには 2 行はさむ:

```lean
haveI : Fintype Q := Fintype.ofFinite Q
rw [Nat.card_eq_fintype_card] at hp
obtain ⟨q, hq⟩ := exists_prime_orderOf_dvd_card (G := Q) p hp
```

`Fintype` が無いまま呼ぶと `rcases` が
`x✝ : ?m is not an inductive datatype` という**原因と無関係な形**で落ちる。

**(d) 仮定に `[Finite (… G …)]` を持つ補題は `(G := G)` を明示しないと
「typeclass instance problem is stuck」になる。**
`obtain ⟨τ, hτ⟩ := exists_orderOf_thetaMul_eq_natCard_quot (A := A) hα hα0 hadj` は
`(A := A)` だけでは `G` がメタ変数のままなので、`[Finite (↥(lowerRamificationGroup B ?G 0) ⧸ …)]`
の解決が止まる。★**インスタンス引数に現れる暗黙引数は全部明示する**（ここでは `(G := G)` も）。

## #116 `def` が包む部分構造は「担い手の型」を明示引数にしないとメタ変数になる（2026-09-07、Y14 / 固定環の塔）

`FixedPoints.subring B ↥H` を薄く包んで

```lean
def fixedRing {B : Type*} [CommRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    (H : Subgroup G) : Subring B := FixedPoints.subring B ↥H     -- ✗
```

と書くと、`fixedRing H` からは **`B` が復元できない**（`H : Subgroup G` は `G` しか決めない）。
その結果、使う側で

```
typeclass instance problem is stuck
  MulSemiringAction G ?m.6
```

が出る。★エラーは `MulSemiringAction` の解決失敗に見えるが、原因は
**`def` の暗黙引数の設計**である。`↥(fixedRing H)` の `0` すら
`OfNat (↥(fixedRing H)) 0 is stuck` になって、症状が 5 箇所に散らばる。

★直し方は **担い手の型を明示引数にする**:

```lean
def fixedRing (B : Type*) [CommRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    (H : Subgroup G) : Subring B := FixedPoints.subring B ↥H     -- ○（以後 `fixedRing B H`）
```

同じ穴は `Submodule` / `Subalgebra` / `Subfield` を部分群やイデアルで添字づけて包むとき
すべてに出る。★**「返り値の型に現れているから推論できる」は嘘** ——
Lean は返り値の型からは暗黙引数を解かない。

**(b) `associated_one_iff_isUnit` の向きは `Associated a 1`。**
`Associated 1 c` が欲しいときは `.symm` が要る:

```lean
exact ⟨0, by rw [pow_zero]; exact (associated_one_iff_isUnit.mpr hcu).symm⟩
```

**(c) `rw [← h]` は等式ゴールの左右両方を書き換える。**
`a * (ρ • x) = a * x` のようなゴールで `rw [← smul_… , ← h]` と繋ぐと、
狙っていない側まで潰れて `unsolved goals` になる。★片側だけ動かしたいときは
`conv_lhs` ではなく **`calc` に開く**のがいちばん短い（Y14 は 1 往復で直った）。

---

## #117 `def` で包んだ部分構造には mathlib のインスタンスが降りてこない（2026-09-07、Y16 / 固定環への作用）

★#116(a) の続き。`fixedRing B H := FixedPoints.subring B ↥H` のように
**mathlib の部分構造を素の `def` で包む**と、

* mathlib が `FixedPoints.subring B ↥H` に付けているインスタンス
  （例: `Mathlib/RingTheory/Invariant/Basic.lean:98` の
  `MulSemiringAction (G ⧸ H) (FixedPoints.subring B H)`、`[H.Normal]` 付き）は
  **`↥(fixedRing B H)` には降りてこない**。インスタンス探索は `instances` 透明度で
  動くので、素の `def` を展開しないからである。

★これは**害と益の両方**である。

* 害: 「mathlib に在るのに使えない」。使いたければ
  `inferInstanceAs (MulSemiringAction (G ⧸ H) ↥(FixedPoints.subring B ↥H))` で
  明示的に橋を架けるか、`fixedRing` を `abbrev` にする。
* 益: **自前のインスタンスを `↥(fixedRing B H)` に立てても mathlib 側と衝突しない**
  （ダイヤモンドが立たない）。Y16 は `MulSemiringAction G ↥(fixedRing B H)` を
  直接構成したが、mathlib の商群版と共存できる。

★★**在庫調査の教訓（10 度目の実証）**: 「`H ⊴ G` なら `G` は `B^H` に作用する」を
語で探しても出ない。**`.cache/mathlib-index.txt` を `FixedPoints\.` で 1 回 grep する**と、
`def` / `instance` / `lemma` が 60 行にまとまって出て、
「在るのは `G ⧸ H` 版であって `G` 版ではない」が **1 往復で確定した**。
★型でも語でもなく、**名前空間で grep する**のがいちばん速いことがある。

★併せて #68 の再確認: そのインスタンスは在るのに `inferInstance` が落ちた。
理由は「無い」ではなく **`Mathlib.RingTheory.Invariant.Basic` を import していない**
（`Algebra.IsInvariant` が `Unknown constant` になるかで 0.01 秒で判別できる）。

★Bash の落とし穴（Lean ではない）: `cat > file << 'EOF'` のヒアドキュメントは
この環境の PreToolUse フックで壊れることがある（`unexpected EOF while looking for
matching '` が出る）。**Lean ファイルは Write ツールで書く**。

## #118 `Nat.cast` 経由の `min` / `1` は `min_eq_left` に食わせる前に `Nat.cast_one` する（2026-09-07、Y15 / Hasse-Arf）

`truncENat_coe : truncENat (k : ℕ∞) r = min (k : ℝ) r` を `k = 1` で使うと、
左辺は `min ((1 : ℕ) : ℝ) (x + 1)`、目標の右辺は `(1 : ℝ)` になる。
この 2 つは `Nat.cast_one` で等しいが**構文的には別**なので、

```lean
exact min_eq_left (by push_cast; linarith)   -- ✗ linarith failed（goal が False になる）
```

と落ちる（`min_eq_left : a ≤ b → min a b = a` の `a` が `((1:ℕ):ℝ)` に固定され、
右辺 `(1:ℝ)` と合わないため、unifier が別の分解を試して壊れる）。**先に潰す**:

```lean
rw [show (1 : ℕ∞) = ((1 : ℕ) : ℕ∞) by rfl, truncENat_coe, Nat.cast_one]
exact min_eq_left (by linarith)              -- ✓ 0.06 秒
```

★同型の罠: `Order.le_of_lt_add_one` の名前つき引数は `a` / `b` ではなく **`x` / `y`**。
`(a := …)` と書くと `Invalid argument name` で落ちる（エラーが候補を出してくれる）。

★`IsDiscreteValuationRing.irreducible_iff_uniformizer` は
`open IsDiscreteValuationRing` が無いと `Unknown identifier`。#68 の変種で、
「mathlib に無い」ではなく「**名前空間を開いていない**」。
`Found/PGC` の DVR まわりは `open IsLocalRing IsDiscreteValuationRing` を既定にする。

## #119 部分環に `Algebra A ↥S` を立てるのは 2 行。★大域インスタンスにせず `letI` で入れる（2026-09-07、Y17 / 固定環の単項生成）

`S : Subring B` が `A` の像を含むとき（`h : ∀ a : A, algebraMap A B a ∈ S`）、

```lean
@[reducible]
def Subring.algebraOfMapsTo (S : Subring B) (h : ∀ a : A, algebraMap A B a ∈ S) : Algebra A ↥S :=
  ((algebraMap A B).codRestrict S h).toAlgebra
```

で立つ。★★**塔は `IsScalarTower.of_algebraMap_eq fun _ => rfl` で出る**（実測、一発）。
`RingHom.codRestrict` は `SubsemiringClass` で一般化されているので `Subring` にそのまま当たる。

★**証明 `h` を引数に取るので大域インスタンスにはできない**（したくもない）。
消費側では **`letI` で証明の中に入れる**。★結論の statement に `Algebra A ↥S` が
現れないなら、それでまったく困らない（`x ∈ Algebra.adjoin A {algebraMap ↥S B ϖ}` のように
`B` の側で述べておけば、`letI` は proof term の中に閉じる）。
★どうしても `↥S` の側で述べたいときは **statement の中に `letI ... ;` を書ける**:

```lean
theorem foo ... :
    letI : Algebra A ↥S := Subring.algebraOfMapsTo S hS
    Algebra.adjoin A ({ϖ} : Set ↥S) = ⊤ := by
  letI : Algebra A ↥S := Subring.algebraOfMapsTo S hS
  ...
```

★これは #117(i)（素の `def` 包みには mathlib のインスタンスが降りてこない）の**逃げ道**でもある。
`fixedRing B H := FixedPoints.subring B ↥H` は型としては `Subring B` なので、
`Subring` に対して立てた `Algebra` / `IsScalarTower` はそのまま当たる。

★クラス型を返す `def` には **`@[reducible]` が要る**（付けないと
`Definition ... of class type must be marked with `@[reducible]`` の警告）。

★併せて: `Module.Finite.of_injective` は `[IsNoetherian R N]` を要求する
**半線形版**（`f : M →ₛₗ[σ] N`）だが、`(IsScalarTower.toAlgHom A C B).toLinearMap` を
渡せば `σ = RingHom.id` で普通に unify する。「中間環が有限生成」はこれ 1 行。

## #120 `if` を割った直後の `Continuous fun r => r` を `simp` は閉じない（2026-09-07、Y18 / 上付き分岐群）

**失敗形**: `truncENat x r = if x = ⊤ then r else min (x.toNat : ℝ) r` の連続性を

```lean
unfold truncENat
by_cases hx : x = ⊤
· simp [hx]        -- ⊢ Continuous fun r => r  が残る（エラー: unsolved goals）
```

`simp` は `if_pos` までは進めるが、**η 展開された恒等写像 `fun r => r` に
`continuous_id` を当ててくれない**（`continuous_id` の左辺は `id` であって `fun r => r`）。

**直し**: `simp only` で `if` だけ潰し、`exact continuous_id` を手で当てる。

```lean
by_cases hx : x = ⊤
· simp only [hx, if_true]; exact continuous_id
· simp only [hx, if_false]; exact continuous_const.min continuous_id
```

★同じ形の和の連続性は `continuous_finsetSum`（`continuous_finset_sum` は 2026-09 に
deprecated）＋ `Continuous.div_const` で 1 行。`φ_G(n) = −1 + (Σ_τ min{i(τ), n+1})/|G|` の
連続性はこれだけで出る（Y18 実測 0.04 秒、一発）。

★併せて: `a ≤ b → 0 < c → a / c ≤ b / c` は名前で引くと
`div_le_div_of_nonneg_right` の引数が `0 ≤ c` だったり `0 < c` だったりして当たらない。
**`gcongr` が一発で閉じる**（`positivity` が `0 < ↑(Nat.card S)` を拾う）。

## #121 `ℕ∞` に `WithTop` の補題を `rw` で当てられない（2026-09-07、Y15b / Hasse-Arf 段 2）

**失敗形**: `htop : ∑ τ : H, ramIndex π' (σ * ↑τ) = (⊤ : ℕ∞)` に対して

```lean
rw [WithTop.sum_eq_top] at htop
-- Did not find an occurrence of the pattern
--   @Eq (WithTop ?m.213) (∑ i ∈ ?m.215, ?m.216 i) ⊤
-- in the target expression
--   @Eq ℕ∞ (∑ τ, ramIndex π' (σ * ↑τ)) ⊤
```

`ℕ∞` は `WithTop ℕ` だが、`rw` の統一は `instances` 透明度なので**畳んだままの
`ℕ∞` を開いてくれない**。

**直し**: `rw` をやめて**項の形**で当て、型を名前つき引数で固定する。

```lean
obtain ⟨τ, -, hτ⟩ := (WithTop.sum_eq_top (M := ℕ)).1 htop
```

★併せて 3 つ:

* `mul_top` は `ℕ∞` には**無い**（`Unknown identifier`）。`ENat.mul_top (h : m ≠ 0)` を使う。
  `WithTop.mul_top` は `[DecidableEq α]` を要求するので `ENat` 版の方が軽い。
* `(Nat.card ↥H : ℕ∞) ≠ 0` は `by exact_mod_cast (Nat.card_pos (α := H)).ne'`。
* `Set G` の積（`(M : Set G) * (H : Set G)`）は **`open Pointwise` が無いと
  `failed to synthesize HMul (Set G) (Set G) ?m`**。Y11 の `herbrand_coe_mul_coe_eq` を
  引くファイルは必ず要る。

## #122 `IsPGroup.of_equiv` は「自分」が第 1 引数（2026-09-07、Y15b）

**失敗形**:

```lean
IsPGroup.of_equiv Subgroup.topEquiv (h1 ▸ isPGroup_lowerRamificationGroup_one …)
-- Type mismatch: has type IsPGroup p ↥⊤ … but is expected to have type ?m ≃* G
```

`of_equiv` / `of_surjective` / `to_quotient` はどれも `variable (hG : IsPGroup p G)` を
**暗黙のセクション変数として第 1 引数に持つ**。ドット記法で書くこと。

```lean
(h1 ▸ isPGroup_lowerRamificationGroup_one (A := A) p huniB hπ'.ne_zero hadj).of_equiv
  Subgroup.topEquiv
```

★`Subgroup.zpowers` の所属を `refine ⟨t, ?_⟩` で開くと、ゴールが
**β 簡約されない `(fun x ↦ ↑σ ^ x) t = ↑g`** の形で出て `rw [← QuotientGroup.mk_zpow]` が
当たらない。`show (QuotientGroup.mk σ : G ⧸ H) ^ t = QuotientGroup.mk g` を 1 行挟む。

## #123 `letI` の下で「引数がラムダだけ」の補題を呼ぶと暗黙型が決まらない（2026-09-07、Y20）

**失敗形**（`haveI` で `SMulCommClass` を入れる 1 行）:

```lean
haveI := smulCommClass_quotient_fixedRing (A := A) (fun σ c => quotientSMul_mk_fixedRing σ c)
-- typeclass instance problem is stuck
--   Subgroup.Normal ?m.172
```

`smulCommClass_quotient_fixedRing {A B G H}` の `B`・`G`・`H` は**引数のラムダからは
決まらない**（ラムダの型は期待型から来るが、`haveI :=` には期待型が無い）。
その結果 `[H.Normal]` の合成が `?m` のまま走って止まる。

**直し方**: 期待型を書くか、**暗黙引数を全部名前で渡す**。後者が短い。

```lean
haveI := smulCommClass_quotient_fixedRing (A := A) (B := B) (G := G) (H := H)
  (fun σ c => quotientSMul_mk_fixedRing σ c)
```

★同じ行を `exact` の引数位置に置くと**期待型があるので何も渡さなくても通る**
（本ファイルの `hq` の 2 度目の出現がそれ）。★「stuck」の語が出たら
**その項に期待型が付いているかどうか**を先に見ること。

★併せて 1 つ（在庫）: **局所環の間の代数写像に沿って剰余体の標数は降りる**。
`IsLocalHom` も剰余体の同型も要らず、`mem_maximalIdeal_of_map_mem`（単元の像は単元、
という易しい向き）と `CharP.charP_iff_prime_eq_zero`（`[Nontrivial R]` が要る）の 2 本で
3 行。`CharP (ResidueField B) p → CharP (ResidueField C) p`。
★分岐（完全分岐か不分岐か）は**1 度も要らない**。

## #124 `AlgEquiv.restrictNormalHom` と `restrictNormalHom_surjective` は `E` と `K₁` の役割が逆（2026-09-07、Y19 / 絶対 Galois 群のフィルトレーション）

```
@AlgEquiv.restrictNormalHom      : ... {K₁} [Algebra F K₁] → (E) [Algebra E K₁] [Normal F E] → Gal(K₁/F) →* Gal(E/F)
@AlgEquiv.restrictNormalHom_surjective : ... {K₁} → (E) [Algebra K₁ E] [Normal F K₁] [Normal F E] → Surjective (restrictNormalHom K₁)
```

**同じ字 `E` が、片方では小さい体（制限先）、もう片方では大きい体（制限元）である。**
「`Gal(K̄/K) ↠ Gal(L/K)` の全射性」を素直に
`AlgEquiv.restrictNormalHom_surjective (F := F) (K₁ := ↥L) (E := E)` の形の
名前付き引数で書こうとして `restrictNormalHom (K₁ := ↥L)` と書くと

```
failed to synthesize instance of type class
  Algebra E ↥L
```

で落ちる（大小が逆なので `Algebra` の向きが合わない）。正しい書き方:

```lean
theorem surjective_restrictNormalHom {F E : Type*} [Field F] [Field E] [Algebra F E]
    [Normal F E] (L : IntermediateField F E) [Normal F L] :
    Function.Surjective (AlgEquiv.restrictNormalHom (F := F) (K₁ := E) (L : Type _)) :=
  AlgEquiv.restrictNormalHom_surjective (F := F) (K₁ := (L : Type _)) E
```

★核が固定化部分群であること（`(restrictNormalHom ↥L).ker = L.fixingSubgroup`）は
`AlgEquiv.restrictNormal_commutes` を 2 回使って 10 行。`ext σ` のあと
`MonoidHom.mem_ker` と `IntermediateField.mem_fixingSubgroup_iff` で開き、
`simp only [AlgEquiv.restrictNormalHom, MonoidHom.mk'_apply, AlgEquiv.one_apply]` で
`restrictNormal` の形に落とすところが要点（`ext x` は `Subtype.ext` まで進むので
`show` で書き直そうとすると `1 x` が残って合わない）。

★併せて 1 つ（同じセッションで踏んだ）: **`∃ P, (IsOpen (P : Set Γ) ∧ P.Normal)` は
型注釈が無いと `P : Γ → Prop` に潰れる**（`(P : Set Γ)` の強制が先に効く）。
出るエラーは `Invalid field 'Normal': ... does not contain 'Function.Normal'` で、
一見なにが起きたか分からない。`∃ P : Subgroup Γ, ...` と書けば直る。

## #125 `SetLike` の台集合の等式は `∈` のゴールに `rw` できない（2026-09-07、Y21 / Hasse-Arf 段 2）

Proposition 6.9 は「`(G_n : Set G) * (H : Set G) = (G_{φ(n)} : Set G)`」という
**台集合の等式**で述べてある（商群を作らずに済ませるため）。`G_{n+1} ≤ H` を入れると

```lean
heq : (↑H : Set G) = ↑(ramificationGroupReal ϖ (herbrandPhi π' H (↑n + 1)))
⊢ x ∈ H
```

になるが、ここで `rw [heq]` は

```
Did not find an occurrence of the pattern ↑H in the target expression  x ∈ H
```

で落ちる。★**ゴールの `x ∈ H` は `SetLike` の membership であって `x ∈ (↑H : Set G)`
ではない**（`Set.mem` に unfold されていない）。直し方は、`Set` の形で `have` を立ててから
`SetLike` の membership に渡すこと:

```lean
have hxmem : x ∈ (H : Set G) := by rw [heq]; exact hc
exact hxH hxmem     -- `hxH : x ∉ H` にそのまま通る（こちらは defeq で通る）
```

★逆向き（`Set` の形のゴールに `SetLike` の仮定を渡す）は `SetLike.mem_coe` で開く。

★併せて 2 つ（同じセッションで踏んだ）:

* **`Nat.dvd_sub'` は消えた。** `Nat.dvd_sub' hb ha : N ∣ b - a` と書くと `Unknown constant`。
  いまは `(Nat.dvd_sub_iff_left hab ha).mpr hb`（`hab : a ≤ b`）。
  ★`exact?` が 0.02 秒で出すので、`Nat` の割り算・引き算まわりは名前を覚えず投げるのが速い。
* **`MulSemiringAction G ↥(fixedRing B H)` には `[H.Normal]` が要る**
  （`fixedRingMulSemiringAction`、Y16）。`G` の作用が固定環に降りるのは `H` が正規のときだけ
  なので当たり前だが、`fixedRing B H` を「ただの部分環」と思って
  `[Fintype H]` だけ書くと `failed to synthesize MulSemiringAction G ↥(fixedRing B H)` が
  **補題の宣言行ではなく本体の `fun r => ...` の位置に出る**ので原因が見えにくい。

## #126 `IntermediateField.fixingSubgroup` の `Normal` は mathlib のインスタンスではない —— 商群を「型に」書けない（2026-09-07、Y19b+c / 不分岐部分の除去）

`L/K` が normal でも `(L.fixingSubgroup : Subgroup Gal(K̄/K)).Normal` は
**インスタンスとして登録されていない**。そのため

```lean
def foo ... : (↥(I ⊔ N') ⧸ L.fixingSubgroup.subgroupOf (I ⊔ N')) ≃* ... := ...
```

のように**宣言の型**に商群が現れると

```
failed to synthesize instance of type class
  Mul (↥(I ⊔ N') ⧸ L.fixingSubgroup.subgroupOf (I ⊔ N'))
failed to synthesize instance of type class
  (L.fixingSubgroup.subgroupOf (I ⊔ N')).Normal
```

が出る。★**証明の中の `haveI` では間に合わない**（型のほうが先に elaborate される）。
直し方は、その `def` より**前**にインスタンスとして置くこと:

```lean
instance normal_fixingSubgroup_of_normal (K : PAdicLocalField p)
    (L : IntermediateField K.carrier K.closure) [Normal K.carrier L] :
    (L.fixingSubgroup : Subgroup K.absGal).Normal := by
  rw [← IntermediateField.restrictNormalHom_ker (K := K.carrier) (L := K.closure) (E := L)]
  infer_instance
```

★根拠は `IntermediateField.restrictNormalHom_ker`（`(restrictNormalHom E).ker = E.fixingSubgroup`）
——「核だから正規」。★`#124` と併せて `restrictNormalHom` 系はこの 2 点で必ず止まる。

★同じ形の在庫: `Subgroup.subgroupOf` の正規性は `Normal` から自動で降りる
（`H.Normal → (H.subgroupOf S).Normal` はインスタンスがある）ので、
足りないのは**いちばん外側の 1 つだけ**である。

## #127 `iff` の前に**明示引数**があると `.mp` が「Unknown constant」になる（2026-09-07、Y19b+c）

```lean
IsLocalRing.residue_eq_zero_iff.mp h
-- Unknown constant `IsLocalRing.residue_eq_zero_iff.mp`
```

`rw [IsLocalRing.residue_eq_zero_iff]` は通るのに `.mp` が通らない。理由:

```
@IsLocalRing.residue_eq_zero_iff : ∀ {R} [CommRing R] [IsLocalRing R] (x : R),
  residue R x = 0 ↔ x ∈ maximalIdeal R
```

★**`(x : R)` が明示引数**なので、この定数の型は `Iff` ではなく `Pi` である。
定数に対する `C.mp` は**名前解決**（`C.mp` という定数を探す）なので落ちる。
（暗黙引数だけなら型は `Iff` に見えるので `.mp` が通る——そこが紛らわしい。）

直し方は `_` を 1 つ入れるだけ:

```lean
(IsLocalRing.residue_eq_zero_iff _).mp h     -- OK
```

★**見分け方**: `#check @foo` して、`↔` の**前**に `(...)` の丸括弧があれば `.mp` は使えない。
★`rw` / `simp` は引数を自分で埋めるので気づかない。**term モードに移した瞬間に出る。**

## #128 「mathlib に無い」の 3 連続誤判定 —— 部分群への作用の制限は**全部ある**（2026-09-07、Y22）

先行ノード（Y21）が「`MulSemiringAction ↥H B` は mathlib にインスタンスが無い」と報告したが、
**誤りだった**。実際は `↥H` への制限に必要なものが 3 つとも在る:

| 要るもの | mathlib の宣言 | 場所 |
|---|---|---|
| `MulSemiringAction ↥H B` | `Subgroup.mulSemiringAction` | `Algebra/Ring/Action/Subobjects.lean:40` |
| `SMulCommClass ↥H A B` | `Subgroup.smulCommClass_left` | `Algebra/Group/Subgroup/Actions.lean:43` |
| `FaithfulSMul ↥H B` | 無名 instance | `Algebra/Group/Subgroup/Actions.lean:59` |

★どれも `inferInstance` で一発で解決する。**自分で `instance` を書いてはならない**
（書くと mathlib のものと 2 本立ちして、`herbrandPhi` のような
`[Fintype ↥H]` を持ち回る定義で**別インスタンスに分岐する**）。

★★**引き方**: `.cache/mathlib-index.txt` を「型」ではなく
**名前空間 + 型クラス名**で grep する（#117(ii)）。1 回で出る:

```
grep -n "FaithfulSMul" .cache/mathlib-index.txt | grep -i "subgroup\|submonoid"
```

★同じ探索で `AddSubgroup.subgroupOf_inertia` も出た:

```
(I.inertia G).subgroupOf H = I.inertia ↥H        -- Algebra/Group/Subgroup/Basic.lean:1077
```

`Ideal.inertia` は `(Submodule.toAddSubgroup I).inertia` の **reducible** な別名なので、
`lowerRamificationGroup B G n := (𝔪^(n+1)).inertia G` と定義してあれば

```lean
theorem lowerRamificationGroup_subtype (H : Subgroup G) (n : ℕ) :
    lowerRamificationGroup B ↥H n = (lowerRamificationGroup B G n).subgroupOf H :=
  (AddSubgroup.subgroupOf_inertia _ H).symm
```

が**定義を展開せずに**通る（`n` について全域、端点の例外なし）。

## #129 `⊤` の場合が既にある補題は、一般の部分群の方が**易しい**ことがある（2026-09-07、Y22）

`herbrandPhiGroup_eq_herbrandPhi_top`（`φ_G = φ_⊤`）は
`Fintype.sum_equiv Subgroup.topEquiv` と `Nat.card_congr` を要していた。
`↥(⊤ : Subgroup G)` と `G` が**別の型**だからである。

ところが一般の `H` については

```lean
theorem herbrandPhi_eq_herbrandPhiGroup_subtype (α : B) (H : Subgroup G) [Fintype ↥H] (n : ℝ) :
    herbrandPhi α H n = herbrandPhiGroup ↥H α n := by
  rw [herbrandPhi_eq_phiOf, herbrandPhiGroup]
  rfl
```

★**`rfl` で閉じる**（0.02 秒）。両辺とも添字型が同じ `↥H` で、
`ramIndex (G := ↥H) α τ` と `ramIndex (G := G) α ↑τ` は定義的に等しいからである。

★教訓: 「`⊤` で苦労した補題だから一般でも苦労する」は**逆**のことがある。
`⊤` の困難は `↥⊤ ≠ G` という**型の食い違い**であって、一般化の困難ではない。
★まず `rfl` を 1 回叩くこと（0.02 秒で終わるか、即座に落ちる）。

## #130 `rw [Nat.cast_zero, Nat.cast_zero]` は 2 つ目で必ず落ちる（2026-09-07、Y22）

`rw` は**その書き換えの全出現**を一度に潰す。`herbrandPhi π H ↑0 = ↑0` のように
`↑0` が左右に 1 つずつある形で `Nat.cast_zero` を 2 回並べると、
1 回目で両方消えるので 2 回目が

```
Tactic `rewrite` failed: Did not find an occurrence of the pattern ↑0
```

で落ちる。★「左右にあるから 2 回」と数えないこと。`simp` に逃げる前に**まず 1 回で試す**。


## #131 `Eq.ge` の `.le` は落ちる —— 等式から割り切りへ渡すのは `▸` の一択（2026-09-07、Y23 / 順分岐商の塔）

`hk : Nat.card ↑N = p ^ k` から `p ∣ Nat.card ↑N` を作りたくて

```lean
dvd_trans (dvd_pow_self p hk0.ne') hk.ge.le
```

と書いたら

```
Invalid field `le`: The environment does not contain `Nat.le.le`,
so it is not possible to project the field `le` from an expression Eq.ge hk
```

`Eq.ge` が返すのは `Nat.le`（構造そのもの）であって `LE.le` の**ドット記法が効く形**
ではない。★等式で型を移すだけなら **`hk ▸ e`** が一番短い:

```lean
(hk ▸ dvd_pow_self p hk0.ne' : p ∣ Nat.card ↑N).trans hdvd
```

★型注釈を付けること（付けないと `▸` の書き換え方向が決まらない）。

## #132 `ENat.one_le_iff_ne_zero` は非推奨。名前空間ごと `Order.` に移った（2026-09-07、Y23）

`ℕ∞` で `0 < x` と `¬ (1 < x)` から `x = 1` を出すとき

```lean
le_antisymm (by exact_mod_cast not_lt.mp hnot) (ENat.one_le_iff_ne_zero.mpr hpos.ne')
```

は**通るが warning が出る**。置き換えは `Order.one_le_iff_ne_zero`。
★★**ドット記法が使えなくなる**点に注意（`x.one_le_iff_ne_zero` →
`Order.one_le_iff_ne_zero x`）。上の形は `.mpr` を付けているだけなので
名前を差し替えるだけで通る。

## #133 `haveI` で入れたインスタンスは**補題を切り出すと消える**（2026-09-07、Y23）

証明を 2 本の宣言に割ったとき、元の証明の中で

```lean
haveI := smulCommClass_quotient_fixedRing (A := A) (B := B) (G := G) (H := H) hq
```

としていたものが、切り出した側にしか残らず、呼び出し側で

```
failed to synthesize instance of type class
  SMulCommClass (G ⧸ lowerRamificationGroup B G 1) A ↥(fixedRing B (lowerRamificationGroup B G 1))
```

になった。★**`letI`/`haveI` は宣言の境界を越えない。**
抽象核と具体層に割るときは、`haveI` を**両方に**書くか、
インスタンス引数として明示的に持ち回ること。
★これは「抽象核をまず切り出す」設計と必ずセットで起きる失敗形である
（1 往復 11 秒で直る。落ちる場所も 1 行で分かる）。

## #134 名前が「無い」4 連発 —— `Subgroup.map_top` / `Subgroup.comap_bot` / `Ideal.pow_le_pow_left` / `IsGalois.toNormal`（2026-09-07、Y19d / 段データの組み立て）

いずれも**綴りが違うだけ**で mathlib にある（#68 の「import していないだけ」とは別種）。

| 書いた名前 | 実際 | 備考 |
|---|---|---|
| `Ideal.pow_le_pow_left` | **`Ideal.pow_right_mono`** | `I ≤ J → I^n ≤ J^n` |
| `Subgroup.comap_bot` | **`MonoidHom.comap_bot`** | `comap f ⊥ = f.ker`。名前空間が `Subgroup` でない |
| `Subgroup.map_top` / `map_top_eq_range` | **無い** | `Subgroup.map H.subtype ⊤ = H` は `simp` も閉じない。`le_antisymm (Subgroup.map_subtype_le _) (fun g hg => ⟨⟨g, hg⟩, Subgroup.mem_top _, rfl⟩)` と手で書く（3 行） |
| `IsGalois.toNormal` | **無い** | `haveI := hgal; infer_instance` で `Normal F E` が出る |

★`exact?` は 0.02〜0.09 秒で返る。**綴りを 3 回試すより 1 回 `exact?` を投げる**方が安い。

## #135 構造体の**引数に依存する項**を `rw [← h]` で書き換えると motive が壊れる（2026-09-07、Y19d）

```lean
structure StageGenerator (K) (N : Subgroup K.absGal) where
  fixingSubgroup_eq : (adjoin K {gen}).fixingSubgroup = N
...
theorem le_stage (g : StageGenerator K N) : N ≤ g.stage v := by
  rw [← g.fixingSubgroup_eq]   -- ★落ちる
```

```
motive is not type correct: fun _a ↦ _a ≤ g.stage v
  g has type StageGenerator K N but is expected to have type StageGenerator K _a
```

★`g` の**型が `N` に依存している**ので、`N` を書き換えると `g` の型も変わる。
**直し方**: ゴールの `N` を書き換えず、**仮定の側**を書き換える。

```lean
  intro σ hσ
  have h1 : σ ∈ (adjoin K {g.gen}).fixingSubgroup := by rw [g.fixingSubgroup_eq]; exact hσ
```

★同じ理由で `calc` を使うと通る場面もある（`rw [← hv]` が `N` を全部書き換えてしまう例）。

## #136 `dite` の条件が `Nonempty` だと `Decidable` が付かない —— `letI := Classical.dec` を**定義の中**に（2026-09-07、Y19d）

```lean
noncomputable def absGalStage (K) (N) (v : ℝ) : Subgroup K.absGal :=
  if h : Nonempty (StageGenerator K N) then h.some.stage v else ⊤
```
→ `failed to synthesize Decidable (Nonempty (StageGenerator K N))`。

★直し方は `letI := Classical.dec (Nonempty (StageGenerator K N))` を**本体の先頭に置く**。
そのうえで、展開する補題は

```lean
  rw [absGalStage]          -- ここまで
  exact dif_pos h           -- ★`rw [absGalStage, dif_pos h]` は通らない
```

と**2 段に割る**（`rw` は `letI` の下の `dite` のインスタンスを合わせられない）。

## #137 Y18 `upperRamification_coe_mul_coe_eq` の右辺は `C` で計算した `G^m` である（2026-09-07、Y19d）

```lean
{ϖ : C} (hϖ : algebraMap C B ϖ = π'') ...
  ↑(upperRamificationGroup G π' m) * ↑H = ↑(upperRamificationGroup G ϖ m)
```

右辺の `ϖ` は **`C = B^H` の元**であり、`upperRamificationGroup G ϖ m` は
**`C` の上で**（`G` の `C` への作用と `C` の付値で）計算した部分群である
——`π'' = algebraMap C B ϖ` を渡すのではない。`π''` は `hfix` にしか出てこない。
★`algebraMap C B ϖ` と書いて型不一致で 1 往復落とした。
★これが「(G/H)^m の引き戻し」の正体で、`H` の商群を作らずに済ませる仕掛けである。

## #138 `Subgroup.card_dvd_of_injective` は**行き先の型で単一化する** —— `Multiplicative X` は `show` で先に据える（2026-09-07、Y24 / Cor 6.13 (iii)）

```lean
-- θ_n : G_n/G_{n+1} ↪ Multiplicative (ResidueField B)
rw [show (…).relIndex (…) = Nat.card (… ⧸ …) from Subgroup.index_eq_card _]
exact Subgroup.card_dvd_of_injective _ (thetaAddQuot_injective (A := A) hα hα0 hadj hi)
```
→ `failed to synthesize Group (ResidueField B)` ＋ `Application type mismatch`。

原因は `Subgroup.card_dvd_of_injective (f : α →* H) : Nat.card α ∣ Nat.card H` の
**`H` が結論の右辺 `Nat.card (ResidueField B)` から先に決まる**こと。
`ResidueField B` は体なので `Group` にはならず、そこで止まる。

★直し方は**行き先を先に据える**:

```lean
show _ ∣ Nat.card (Multiplicative (ResidueField B))
exact Subgroup.card_dvd_of_injective _ (thetaAddQuot_injective (A := A) hα hα0 hadj hi)
```

`Multiplicative X` は型シノニムなので `Nat.card (Multiplicative X) = Nat.card X` は
定義的に等しく、`show` 1 行で通る（`rfl` 補題を探しに行かなくてよい）。
★同じ形は `Additive` 側でも起きる。

## #139 `Subgroup.relIndex` は `Subgroup.index_eq_card` で商の `Nat.card` になる（2026-09-07、Y24）

`H.relIndex K` の定義は `(H.subgroupOf K).index` なので、

```lean
show (lowerRamificationGroup B G 1).relIndex (lowerRamificationGroup B G 0)
  = Nat.card (lowerRamificationGroup B G 0 ⧸
      (lowerRamificationGroup B G 1).subgroupOf (lowerRamificationGroup B G 0)) from
  Subgroup.index_eq_card _
```

が**そのまま通る**（`Subgroup.relIndex` を `unfold` しなくてよい）。
★この木の `thetaMulQuot` / `thetaAddQuot`（Prop 6.2）の定義域はまさにこの商なので、
**分岐群の商の位数は `relIndex` で書くのがいちばん安い**。
望遠鏡積 `[G : F n] = ∏_{i<n} [F i : F (i+1)]` も
`Subgroup.relIndex_mul_index : H ≤ K → H.relIndex K * K.index = H.index` の帰納 4 行で出る
（★`mul_comm` を 1 つ挟む必要がある —— `Finset.prod_range_succ` の向きが逆）。

## #140 `Polynomial.Monic` の項に**ドット記法は使えない** —— `Eq` の名前空間へ落ちる（2026-09-07、Y25 / Λ6 §4-a）

`Polynomial.Monic p` は `p.leadingCoeff = 1` の `def` である。したがって

```lean
rw [(hmonic.map φ).natDegree_eq_of_map] at hk   -- ✗
```

は `Polynomial.Monic.natDegree_eq_of_map` ではなく **`Eq.natDegree_eq_of_map`** を探しに行き、

```
Invalid field `natDegree_eq_of_map`: The environment does not contain `Eq.natDegree_eq_of_map`,
so it is not possible to project the field `natDegree_eq_of_map` from an expression
  Monic.map φ hmonic
of type `(map φ g).leadingCoeff = 1`
```

という**行き先の名前空間が違う**エラーになる（「そんな補題は無い」ではない）。
★直し方は「元の `Monic` 項に補題を当てる」:

```lean
have hdegmap : (g.map φ).natDegree = g.natDegree := hmonic.natDegree_map φ   -- ○
```

`Polynomial.Monic.natDegree_map (hmo : P.Monic) (f : R →+* S) : (P.map f).natDegree = P.natDegree`。
★同じ理由で `(hmonic.map φ).leadingCoeff` は**書ける**（`Monic.leadingCoeff` は
`Monic` の項を取る補題で、ドット記法が `Eq` に落ちても引数として渡るため）。
★★一般則: **`def` で `Prop` に展開される述語（`Monic` / `IsUnit` ではない方）**の項に
ドット記法を使うときは、展開先の名前空間（ここでは `Eq`）が優先されると思ってよい。

## #141 `lean/ABC3/Found.lean` は **CRLF** —— Python の文字列置換は `'rb'` で読むこと（2026-09-07、Y25）

`Found.lean` に import を 1 行足すとき、

```python
s = io.open(p, encoding='utf-8').read()      # ← 既定は universal newlines
assert s.count('import ABC3.Found.PGC.Foo\n') == 1   # ✗ 0 になることがある
```

は `newline=''` を付けると `\r\n` のまま入るので `\n` で数えると **0 件**になり、
逆に `newline` 既定だと読めても**書き戻しで CRLF が LF に潰れて 1786 行全部が diff に出る**。
★安全なのは **binary で読んで binary で書く**:

```python
raw = io.open(p, 'rb').read()
eol = b'\r\n' if raw.count(b'\r\n') > 0 else b'\n'
old = b'import ABC3.Found.PGC.Bar' + eol
raw = raw.replace(old, old + b'import ABC3.Found.PGC.Foo' + eol)
io.open(p, 'wb').write(raw)
```

★★`git diff --stat lean/ABC3/Found.lean` が **1 行**であることを必ず確認する
（並行セッションが足した行が混ざるので `git diff` の中身も見ること）。
★`cat >> f << 'EOF'` / `python - <<'PYEOF'` はどちらも PreToolUse フックに
潰される（#117(iv)）。**`.py` を Write して呼ぶ**のがいちばん速い。

## #142 `set` は**他の仮説の型に現れる項**を抽象化すると、その仮説を `τ✝` に化けさせる（2026-09-07、Λ6 §4-b）

```lean
(τ : ↥(IntermediateField.adjoin F' ({φ x} : Set M)) ≃ₐ[F'] ↥(IntermediateField.adjoin F' ({φ x} : Set M)))
...
set B := IntermediateField.adjoin F' ({φ x} : Set M) with hB   -- ✗
```

`set` は**ゴールと文脈の両方**を書き換えるので、`τ` の型も `↥B ≃ₐ[F'] ↥B` に変わる。
このとき Lean は `τ` を**新しい局所変数に取り替え**、元の `τ` は `τ†` として残る。
結果、ゴールに出てくるのは `τ†` の方で、`exact hsx`（`hsx` は新しい `τ` の言明）が

```
Type mismatch: has type ↑(τ† ⟨φ y, ⋯⟩) = φ (s y) but is expected to have type ↑(τ ⟨φ y, ⋯⟩) = φ (s y)
```

で落ちる。★**束縛変数の型に現れる項に `set` を使わない。** 長くて読みにくくても
`have hkey : <長い式> = <長い式> := by ...` と**その場に書き下す**方が速い（実測 1.72 秒 → 0.26 秒）。

## #143 `∃!` に `refine ⟨_, ?_, ?_⟩` するとβ簡約されず、`rw` が「パターンが見つからない」と言う（2026-09-07、Λ6 §4-b）

`∃! ρ, P ρ` は `∃ ρ, P ρ ∧ ∀ y, P y → y = ρ` なので、`refine ⟨a, ?_, ?_⟩` の第 1 ゴールは

```
⊢ (fun ρ ↦ closureCompletionCoe K ↑(ρ ⟨x, hmem⟩) = ↑(τ ⟨…⟩)) a
```

という**β簡約されていない**形で出る。ここで `rw [coe_algEquivRestrictSelf]` は
`Did not find an occurrence of the pattern` で落ちる（パターンはラムダの中にある）。

★直し方は 2 つ。**`have hmain : <明示的な型> := by …` を先に作って `refine ⟨a, hmain, ?_⟩`**
（これがいちばん安全。型を書くので以後 `rw` が効く）か、`dsimp only` でβ簡約する。
★一意性側（`intro ρ hρ`）は `intro` がβ簡約するので**そのままで通る**（非対称なので注意）。

## #144 複数行にまたがる `calc` の第 1 項が関数適用だと、パーサが途中で切る（2026-09-07、Λ6 §4-b）

```lean
calc closureCompletionCoe K
    ((… : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) : K.closure)
    = closureCompletionCoe K (…) := by …
  _ = … := …
```

は `invalid 'calc' step, left-hand side is … but previous right-hand side is closureCompletionCoe K`
（＝**第 1 項が `closureCompletionCoe K` だけだと解釈された**）で落ち、続けて
`unexpected token '_'` になる。★**長い `calc` は書かず、`have h1 : … := …` を並べて
最後に `rw […]; exact …` で閉じる**方が安い（実測: `calc` 4 段を `have` 1 本 + `rw` 1 本に潰せた）。

## #145 `adjoinIntegersIncl` の係数を 2 層いっぺんに `rfl` で潰すと kernel が落ちる（2026-09-07、Y19e）

`#full-check`

```lean
-- ★NG: (kernel) deterministic timeout
theorem coeC (hle) (z : adjoinIntegers K x) :
    (((adjoinIntegersIncl K hle z : adjoinIntegers K x') : ↥K.carrier⟮x'⟯) : K.closure)
      = ((z : ↥K.carrier⟮x⟯) : K.closure) := rfl
-- ★NG も同じ: ((adjoinIntegersIncl K hle z).1.1 : K.closure) = (z.1.1 : K.closure) := rfl
```

★**1 層なら通る**（実測 0.1 秒台）:

```lean
theorem val_adjoinIntegersIncl (hle) (z : adjoinIntegers K x) :
    ((adjoinIntegersIncl K hle z).1 : ↥K.carrier⟮x'⟯) = ⟨(z.1 : K.closure), hle z.1.2⟩ := rfl
```

2 層目は **`congrArg Subtype.val` で上げる**（`rfl` で書かない）。
★`TotallyRamified.lean` の `norm_mk_of_le` が同じ回避を先にやっている（#59 の親戚）。

## #146 `ABC3.Found.PGC.ker_restrictNormalHom_eq_fixingSubgroup` は同名が 2 つある（2026-09-07、Y19e）

`AbsGalRamificationFiltration.lean:492`（`(L : IntermediateField F E)` が最後の明示引数）と
`LubinTateClosure.lean:127`（`(E : IntermediateField F Ω)` が明示引数）。
**片方しか import されていないうちは動くが、両方が import 圏に入ると曖昧になって落ちる**
（`AbelianSplitUnramified` を足した瞬間に `Application type mismatch` が出た）。

★直し方: mathlib の **`IntermediateField.restrictNormalHom_ker (K := …) (L := …) (E := …)`** を使う。
★同じ罠は他の重複名にもある。`grep -c "^theorem <名前>" ` ではなく
`grep -n "PGC\.<名前>" .cache/decl-index.txt` で**何個あるか**を先に見ること。

## #147 この mathlib の `Subgroup.mul_normal` は `↑(H ⊔ N) = ↑H * ↑N` の向き（2026-09-07、Y19e）

`rw [Subgroup.mul_normal]` は `↑(H ⊔ N)` を探し、`rw [← Subgroup.mul_normal]` が `↑H * ↑N` を探す。
★教科書の記憶（`↑H * ↑N = ↑(H ⊔ N)`）と逆なので、**1 往復無駄にしやすい**。
`Did not find an occurrence of the pattern ↑(?H ⊔ ?N)` が出たら向きを疑うこと。

## #148 `∀ U ∈ s, P U` の穴に、`{U}` を暗黙にした補題は嵌まらない（2026-09-07、Y19f）

近傍基の仮説

```lean
(hbasis : ∀ U ∈ nhds (1 : Γ), ∃ N ∈ F.base, (N : Set Γ) ⊆ U)
```

に、こう書いた補題を渡すと落ちる:

```lean
theorem foo (K) (hcompat) {U : Set K.absGal} (hU : U ∈ nhds 1) : ∃ N ∈ …, ↑N ⊆ U
```

```
Application type mismatch: the argument `foo K hcompat` has type
  ?m ∈ nhds 1 → ∃ N ∈ …, ↑N ⊆ ?m
but is expected to have type
  ∀ U ∈ nhds 1, ∃ N ∈ …, ↑N ⊆ U
```

★`U` が暗黙だと**先頭の `∀ U` が消える**ので、`∀ U ∈ …` の形と合わない。
直し方は 2 つ: 補題側の `U` を**明示引数にする**（推奨。1 文字の差）か、
呼び出し側で `fun U hU => foo K hcompat hU` と η 展開する。
★「近傍基」「開被覆」など `∀ U ∈ …` を仮説に取る抽象核を書くときに必ず当たる。

## #149 `π` は `hπmax` の型に現れるので、`rw [← h]`（`h : … = π`）は必ず `motive is not type correct` になる（2026-09-07、Λ6 §4-c）

Lubin-Tate の木では素元 `π` が

```lean
{π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
```

の形で**すべての宣言の添字**になっている。したがって `h : a * u = π` を得ても

```lean
rw [← h]   -- ✗ motive is not type correct:
           --   fun _a => … = _a * ↑u⁻¹ の中で hπmax の型が壊れる
```

は通らない（`iteratedLubinTatePsi hq hπmax …` が `hπmax : … = span {_a}` を要求してしまう）。

★直し方は「**書き換えずに済ませる**」——`Units` の移項補題をそのまま使う:

```lean
exact ⟨u⁻¹, (Units.eq_mul_inv_iff_mul_eq u).mpr h⟩
```

（`Units.eq_mul_inv_iff_mul_eq (c : αˣ) {a b} : a = b * ↑c⁻¹ ↔ a * ↑c = b`。
★`c` が**明示**・`a b` が暗黙なので `(… u).mpr h` の順で書く。`.mpr` を先に書くと
`Invalid field 'mpr'` になる。）
★同じ形は `π` を右辺に持つ等式すべてに出る。`conv` や `nth_rewrite` で逃げるより、
移項補題を探す方が速い（実測 1 往復）。


## #150 抽象核の「不変性」仮説を `∀ (c : M) {x : A}, x ∈ s → c • x ∈ s` と書くと、代入側の `fun σ hx => …` が暗黙引数に食われる（2026-09-07、Y25）

不変部分環へ作用を制限する抽象核

```lean
def mulSemiringActionOfSubringClass … (s : S)
    (hs : ∀ (c : M) {x : A}, x ∈ s → c • x ∈ s) : MulSemiringAction M ↥s
```

に対して、具体層で

```lean
mulSemiringActionOfSubringClass (absClosureInt K) (fun σ hx => smul_mem_absClosureInt K σ hx)
```

と書くと、`hx` が**暗黙の `{x : A}`** に束縛されて

```
Application type mismatch: The argument hx has type `K.closure` of sort `Type`
but is expected to have type `?m ∈ s` of sort `Prop`
```

となる。★悪いのはさらに先で、instance の定義が壊れたまま下流へ伝播し、
`↑(σ • w) = σ ↑w` が `rfl` で閉じない（"Not a definitional equality"）という
**別の顔のエラーが 4、5 個同時に出る**。原因は 1 行目だけである。

★直し方: **抽象核側で `x` を明示にする**。

```lean
    (hs : ∀ (c : M) (x : A), x ∈ s → c • x ∈ s)
  smul c x := ⟨c • (x : A), hs c x x.2⟩        -- 呼び出し側は `fun σ _ hx => …`
```

★教訓：**抽象核の仮説に暗黙引数を置かない**。
抽象核は「代入するだけ」にするのが目的なので、仮説の引数は全部明示の方が安い。
（#148 と同じ味である——「穴を暗黙にした抽象核は嵌まらない」。）

## #151 section の `variable` 仮説（`hπne0`・`hπmax`）は「文に現れない」と黙って落ちる（2026-09-07、Y26 / Lemma 4.3(ii)）

**症状**: `variable (hπne0 : π ≠ 0)` を section に置いてある。
文（statement）には `π` しか出てこないが、証明の中で
`LubinTateAction_comp hq hπmax hπne0 …` と使う。すると

```
error: Unknown identifier `hπne0`
```

が**証明の中で**出る。Lean 4 の自動 `variable` 取り込みは
**「statement に現れる変数」だけ**を入れるので、証明でしか使わない仮説は入らない。
★同じことが「`π` は現れるが `hπmax` は現れない」補題（例:
`𝔭^m = Ideal.span {π^m}`）でも起きる ——
`hπmax` が落ち、続けて `Application type mismatch: hπmax … expected ℕ` という
**呼び出し側の別の顔のエラー**になる（引数がずれるため）。

★直し方: `include hπne0 in` を `theorem` の直前に置く。
`omit` と併用するときは**この順で並べる**（docstring より前。#107）:

```lean
omit [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])] in
include hπne0 in
/-- docstring -/
theorem foo … := …
```

★もう 1 つの手（今回はこちらを採った箇所もある）: **その補題だけ section の外に出して
引数を明示で書く**。`𝔭^m = (π^m)` のような「π と hπmax しか要らない」補題は、
section 変数の取り込み規則と戦うより `(K) {π} (hπmax) (m)` と書いた方が短い。

## #152 `K.carrier⟮x⟯` は `open IntermediateField` が要る（pretty-printer の出力は貼り戻せない）（2026-09-07、Y26）

**症状**: `#check` の出力が `x ∈ K.carrier⟮x⟯` と表示するのでそのまま書くと

```
error: expected token
```

が `⟮` の位置に出る。この木の定型の前置き
`open ABC3.Skeleton.PGC` / `open scoped NormedField Valued Classical` では
**`⟮⟯` の記法が入らない**（`IntermediateField` 名前空間に scoped されている）。

```lean
open ABC3.Skeleton.PGC IntermediateField in   -- ← これなら通る（実測 0.02 秒）
example … : x ∈ K.carrier⟮x⟯ ↔ True := by simp [IntermediateField.mem_adjoin_simple_self]
```

★`Found/PGC/AdjoinIntegers.lean` 系はどこでも
`IntermediateField.adjoin K.carrier ({x} : Set K.closure)` と**展開して**書いている。
記法を使わない方が、`variable` を跨ぐときの表示ゆれも起きない。

## #153 「`x` を含む有限次 Galois 中間体」は `FiniteGaloisIntermediateField.adjoin` が最短（2026-09-07、Y27 / Thm 6.15 B2+B3）

**用途**: 副有限な `Γ_K` で「`τ = 1` を示す」とき、`τ x = x` を各 `x` について示す。
そのために「`x` を含む**有限次 Galois**な中間体 `M`」と「`σ|_M` の位数」が要る。

`normalClosure K.carrier (adjoin K.carrier {x}) K.closure` と書くと
`FiniteDimensional` / `Normal` の instance を手で並べることになる。
mathlib には**束ねた型**がある（`FieldTheory/Galois/GaloisClosure.lean`）:

```lean
let M : IntermediateField k E :=
  (FiniteGaloisIntermediateField.adjoin k ({x} : Set E)).toIntermediateField
have hxM : x ∈ M := FiniteGaloisIntermediateField.subset_adjoin k ({x} : Set E) rfl
haveI : FiniteDimensional k M := inferInstance   -- ★4 つとも inferInstance で出る
haveI : Normal k M := inferInstance
haveI : Finite ((M : Type _) ≃ₐ[k] (M : Type _)) := inferInstance
```

★`[IsGalois k E]` と `[Finite s]` が要る（`{x}` は自動）。実測 0.30 秒で 4 つとも通った。
★これで `e := orderOf (AlgEquiv.restrictNormalHom (M : Type _) σ)` が `≠ 0`
（`(orderOf_pos _).ne'`）になり、`σ ^ e ∈ M.fixingSubgroup` は
`IntermediateField.restrictNormalHom_ker M` + `MonoidHom.mem_ker` + `map_pow` で出る。

★**開部分群を `Γ_K` 側で作るときは `Subgroup.comap` を使う**（`IntermediateField.map` を
経由しない）。`(unramLevel K e).fixingSubgroup |>.comap (AlgEquiv.restrictNormalHom …)` は
`InfiniteGalois.restrictNormalHom_continuous` で開性がそのまま引き戻せるので、
**中間体の 2 層をまたぐ `rfl`（#59）に一度も触らずに済む**。

## #154 `omit [X] in` は docstring の**前**に置く（#107 の仲間）（2026-09-07、B4 / Prop 6.14 n=1）

`section` の `variable` に `[Fintype G]` があり、ある定理だけそれを使わないとき、
linter が `omit [Fintype G] in` を勧める。だが**docstring の後ろに置くと構文エラー**になる:

```lean
/-- ... -/
omit [Fintype G] in                 -- ✗ unexpected token 'omit'; expected 'lemma'
theorem foo ...
```

正しくは `open … in`（#107）と同じで **docstring の前**:

```lean
omit [Fintype G] in
/-- ... -/
theorem foo ...
```

★同じ波で当たった**名前の改称 2 件**（どちらも `Unknown constant`）:

| 書いたもの | 今の名前 |
|---|---|
| `Nat.pos_pow_of_pos` | `Nat.pow_pos`（引数は底の正値性 1 つだけ） |
| `Nat.Ico_succ_right` | 無い。`Finset.Ico_succ_right_eq_Icc` はあるが `Order.succ` 経由で使いにくい |

★`Finset.Icc 1 (N-1) = Finset.Ico 1 N`（`N ≥ 1`）は**補題を探さず**
`ext n; simp only [Finset.mem_Icc, Finset.mem_Ico]; omega` が最短（`N ≥ 1` は文脈にあればよい）。

★**`omega` は「切り詰め引き算 × 変数」を原子化しない**。
`((r+1)^(M+2) - (r+1)^(M+1)) * K = (r*(r+1)^(M+1)) * K` は omega が落ちる
（`* K` ごと 1 つの原子にしてしまう）。
**引き算だけの等式を先に `have` で出して `rw` してから**掛け算に触ること。

## #154 `AlgEquiv.restrictNormalHom` を含む目標に `exact` を当てると whnf / isDefEq が 200000 heartbeats で落ちる —— `set φ := …` で局所定数にする（2026-09-07、Y27 / Thm 6.15 B2+B3）

**症状**（実測 2 往復ぶん無駄にした）:

```lean
haveI := normal_of_isUnramifiedAdjoin K x hx
haveI := isCyclic_gal_of_isUnramifiedAdjoin K x hx
rw [← IntermediateField.restrictNormalHom_ker, MonoidHom.mem_ker, map_mul, map_mul,
  map_inv, map_inv]
exact commutatorElement_eq_one_of_isCyclic _ _
--    ^ error: (deterministic) timeout at `whnf` / `isDefEq`
```

`(G := …)` で型を明示しても**変わらない**（別の顔で `isDefEq` に出るだけ）。

**直し**: `AlgEquiv.restrictNormalHom …` を `set` で**局所定数**にしてから `rw`／`exact` する。
`Found/PGC/AbelianFrobeniusSplit.lean:186` が同じ形をこう書いていた:

```lean
set Kx := IntermediateField.adjoin K.carrier ({x} : Set K.closure) with hKx
set φ  := AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) (Kx : Type _) with hφ
have hker : φ.ker = Kx.fixingSubgroup := IntermediateField.restrictNormalHom_ker Kx
rw [← hker, MonoidHom.mem_ker]
simp only [map_mul, map_inv]
exact commutator_eq_one_of_mul_comm (mul_comm_of_forall_mem_zpowers hg0 (φ a) (φ b))
```

★**単元化の対象を `φ a` / `φ b` と書き下す**のも効いている（`_ _` にすると
`restrictNormalHom` の codomain を再構成しにいく）。

★ついでに:**`IsCyclic G` のインスタンスを使うより、生成元 `g0` を
`IsCyclic.exists_generator` で出して `mul_comm_of_forall_mem_zpowers hg0` を使う方が安い。**
`IsCyclic.commGroup` を `haveI` すると `Mul` が二重になって
`Type mismatch … instHMul … this.toCommMonoid.toCommSemigroup.toCommMagma.toMul` が出る。

★**在庫に同じ形が既に在るか**を先に見ること。`mul_comm_of_forall_mem_zpowers` /
`commutator_eq_one_of_mul_comm` は `AbelianFrobeniusSplit.lean` の §0 に純群論として
切り出されていた（`grep -n "^theorem" …` で 10 秒）。

## #155 `⁅a, b⁆`（**元**の交換子）は `Bracket G G` のインスタンスが立っておらず使えない（2026-09-07、Y27）

`commutator G`・`Subgroup.commutator_le`・`commutator_def` は**使える**のに、
元の交換子だけは

```
failed to synthesize instance of type class
  Bracket G G
```

になる（`#synth Bracket (Multiplicative ℤ) (Multiplicative ℤ)` でも同じ。
`#check @commutatorElement` は `{G} → [Group G] → Bracket G G` を返すので**定義は在る**が
インスタンスとして登録されていない）。

**直し**: `a * b * a⁻¹ * b⁻¹` と書き下す。`Subgroup.commutator_le` の右辺
（`∀ g₁ ∈ H₁, ∀ g₂ ∈ H₂, ⁅g₁, g₂⁆ ∈ H₃`）には
`show a * b * a⁻¹ * b⁻¹ ∈ H₃` で入れる（定義展開で通る、実測 0 秒）。

## #156 `Subsingleton (↥(⊥ : IntermediateField F E) ≃ₐ[F] ↥⊥)` は**推論されない**（`Subsingleton (F ≃ₐ[F] F)` は推論される）（2026-09-07、B1 / Thm 6.15）

「`A = ⊥` なら `Gal(A/k)` は自明」を退化の自己検査に使いたいとき、素直に

```lean
example {F E : Type*} [Field F] [Field E] [Algebra F E] :
    Subsingleton ((⊥ : IntermediateField F E) ≃ₐ[F] (⊥ : IntermediateField F E)) := by
  infer_instance
```

は落ちる（`failed to synthesize … Subsingleton Gal(↥⊥/F)`、実測 0.06 秒）。
★**同じセッションで `Subsingleton (F ≃ₐ[F] F)` は `infer_instance` が一発で通る**ので、
「`⊥ ≃ₐ F` なのだから同じだろう」という予想は外れる（`IntermediateField.botEquiv` を
自分で挟まないと繋がらない）。

**直し**（`botEquiv` を挟むより短い。実測 0.14 秒・一発）:

```lean
theorem algEquiv_eq_one_of_eq_bot {A : IntermediateField k E} (h : A = ⊥) (σ : A ≃ₐ[k] A) :
    σ = 1 := by
  subst h
  refine AlgEquiv.ext fun x => ?_
  obtain ⟨a, ha⟩ := (IntermediateField.mem_bot).mp x.2
  have hx : x = algebraMap k (⊥ : IntermediateField k E) a := by
    apply Subtype.ext; simpa using ha.symm
  rw [hx]
  simp
```

★**#155 の続き**（元の交換子 `⁅a,b⁆` が書けない件）に効く定型が 2 つある。

1. **自分で `⁅⁆` を書かず、補題に産ませる**と `commutatorElement_*` 系は使える。
   `rw [MonoidHom.mem_ker, map_commutatorElement, commutatorElement_eq_one_iff_mul_comm]`
   は通る（`⁅⁆` は `map_commutatorElement` が作る）。
   同様に `Subgroup.commutator_le.mp h a ha b hb` / `Subgroup.commutator_mem_commutator ha hb`
   は `⁅a,b⁆ ∈ H` を**返して**くれるので、こちらから書かなくてよい。
2. **`x * y * x⁻¹ * y⁻¹ = 1` から `x * y = y * x` へは
   `rwa [mul_inv_eq_one, mul_inv_eq_iff_eq_mul]`**（0.02 秒）。
   `commutatorElement_eq_one_iff_mul_comm` を当てようとすると
   `typeclass instance problem is stuck: Group ?m` になる
   （`⁅⁆` の形に戻せないので暗黙の `G` が決まらない）。

## #157 `IntermediateField.lift_bot` / `lift_top` は明示引数の形が読めない —— `simp` を使う（2026-09-07、B5 / Thm 6.15）

`InfiniteGalois.restrict_fixedField` で `E ⊓ L = lift (fixedField …)` に持ち込んだあと、
`fixedField … = ⊥` を代入して残るのは `IntermediateField.lift ⊥ = ⊥`。
索引の statement は

```
theorem IntermediateField.lift_bot (K : IntermediateField F E) : lift (F := K) ⊥ = ⊥
```

なので「中間体を渡せばよい」と読めるが、**実際に渡すと外側の体だと解釈される**:

```lean
exact IntermediateField.lift_bot (unramifiedClosure K)
-- typeclass instance problem is stuck: Algebra ↥(unramifiedClosure K) ?m
exact IntermediateField.lift_bot _
-- 同じく stuck（?m が 2 つ残る）
```

`lean_check` で型を見ると `lift_bot ↥M : ∀ (K : IntermediateField ↥M ?m), lift ⊥ = ⊥` で、
第 1 明示引数は**下の層の体**（`↥M`）だった。

**直し**: `simp` で閉じる（0.08 秒・一発）。`lift_top` も同様。
★**教訓**: `.cache/mathlib-index.txt` の statement は section variable 名をそのまま出すので、
`(K : IntermediateField F E)` の `K` が「自分の書いている `K`」と一致するとは限らない。
迷ったら `lean_check` に 2 行の `example` を投げて**型を見る**（0.08 秒）。

## #158 「もう在る補題」を最速で見つける手 —— **書いて `already been declared` を出させる**（2026-09-07、B5 / Thm 6.15）

B5 は `fixingSubgroup_sup`（`(E ⊔ F).fixingSubgroup = E.fixingSubgroup ⊓ F.fixingSubgroup`）と
`isTotallyRamifiedAdjoin_of_inf_unramifiedClosure_eq_bot` を「無いだろう」と思って書き始めたが、
**2 本とも木に在った**（前者は `Found/PGC/CompositumSurjection.lean` と
`Found/PGC/AbelianDecomposition.lean` に**2 つ**、後者は `Found/PGC/AbelianSplitUnramified.lean`）。

見つけ方は grep ではなく、**同じ名前で書いて `leanfile.mjs` に投げた**だけである:

```
error: `ABC3.Found.PGC.fixingSubgroup_sup` has already been declared
```

★**名前を当てに行くのは無料**（証明を書く前に `theorem 名前 : 型 := sorry` でもよい）。
★しかも在庫のほうが**仮定が少ない**ことがある ——
`isTotallyRamifiedAdjoin_of_inf_unramifiedClosure_eq_bot` は在庫版が
`[FiniteDimensional]` も `[Normal]` も要らない形だった。
★★**名前が思い付く補題は、まず名前で衝突させて確かめること。**

## #159 `adjoinIntegers K x` のノルムを `spectralNorm` の補題と繋ぐには `simpa` ではなく `show`（2026-09-07、B4 段 1 / Prop 6.14）

`y : adjoinIntegers K x` の `‖y‖` は `spectralNorm K.carrier K.closure ↑↑y` と**定義的に等しい**が、
`spectralNorm_…` 系の補題を `simpa using` で当てると**必ず落ちる**:

```
Type mismatch: After simplification, term `this` has type
  spectralNorm K.carrier K.closure ↑↑(…) = ‖↑π‖ ^ (↑e)⁻¹
but is expected to have type
  ‖↑↑(…)‖ = ‖↑π‖ ^ (↑e)⁻¹
```

`simp` が `‖·‖` を**畳まずに** `↑↑` の側だけ動かすので、両辺が別の顔になる。

**直し**: `show` で `spectralNorm` の顔に**先に**書き換えてから `exact`（`simpa` を使わない）:

```lean
show spectralNorm K.carrier K.closure
  (↑(↑(EXPR : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)) = _
exact spectralNorm_of_mem_iteratedLubinTatePsiTorsionPoints K … 
```

★同じ形は `‖α‖ < 1` にも要る（`spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints`）。
★`Found/PGC/LubinTateActionInjective.lean` が既にこの定型を使っていた（`show spectralNorm … < 1`）。

## #160 部分型の `.2` に型注釈を付けると壊れる（`(ϖ.2 : ‖ϖ‖ ≤ 1)` は不可）（2026-09-07、同上）

`adjoinIntegers K x` は `{y | ‖y‖ ≤ 1}` なので `ϖ.2 : ‖ϖ‖ ≤ 1` が**通る文脈もある**が、
`rcases lt_or_eq_of_le (ϖ.2 : ‖ϖ‖ ≤ 1)` のように**型注釈の中**に置くと

```
Application type mismatch: `ϖ.property` has type `↑ϖ ∈ adjoinIntegers K x`
but is expected to have type `?m ≤ ?m`
```

**直し**: `have hϖle : ‖ϖ‖ ≤ 1 := ϖ.2` と**独立した `have` にする**（期待型から unfold が走る）。
★`fun a => a.2` を `(hle : ∀ a : B, N a ≤ 1)` に渡すのは**通る**（期待型が先に決まるため）。

## #161 名前が無い 2 つ ——`isUnit_of_mul_eq_one` / `IsLocalRing.mem_maximalIdeal`（2026-09-07、同上）

* `isUnit_of_mul_eq_one` は **無い**。`isUnit_iff_exists_inv.mpr ⟨c, h⟩` を使う。
* `IsLocalRing.mem_maximalIdeal` も **無い**。`maximalIdeal R` は `nonunits R` と定義的に等しいので
  **`mem_nonunits_iff.mp` / `.mpr` が両向きに使える**
  （`h : π ∈ maximalIdeal 𝒪` から `mem_nonunits_iff.mp h : ¬IsUnit π`、逆も同様）。
* `a ^ k = a ^ n → k = n`（実数、`0 < a < 1`）は `exact?` が**見つけられない**。
  `(pow_right_strictAnti₀ h0 h1).injective` を明示的に組む（`exact?` は `StrictAnti` の形なら即答する）。


## #162 `Subgroup.comap_bot` は無い / `MonoidHom.ker_eq_bot_iff` は**引数を明示**しないと `.2` が使えない（2026-09-07、B5 の穴 1 / Cor 6.13(i)）

`comap f ⊥ = f.ker` の名前は `Subgroup.comap_bot` **ではなく** `MonoidHom.comap_bot`
（`Algebra/Group/Subgroup/Ker.lean`）。相方は `QuotientGroup.ker_mk' : (mk' N).ker = N`。
この 2 本で `Subgroup.comap (QuotientGroup.mk' H) ⊥ = H` が `rw` 一発で閉じる。

さらに `MonoidHom.ker_eq_bot_iff` は **`(f : G →* N)` を明示引数に取る**ので

```lean
MonoidHom.ker_eq_bot_iff.2 e.injective        -- ✗
```

は

```
Invalid projection: Projections cannot be used on functions, and
  MonoidHom.ker_eq_bot_iff has function type ∀ (f : ?m →* ?m), f.ker = ⊥ ↔ Function.Injective ⇑f
```

で落ちる。★**直し**は引数を書く: `(MonoidHom.ker_eq_bot_iff (e : G →* G')).2 e.injective`。
★同じ顔のエラー（`Invalid projection`）は「`.2` の前の項がまだ関数」のときに必ず出るので、
**`iff` の補題に明示引数が残っていないか**を最初に疑うこと。

## #163 `Fintype.sum_equiv e.toEquiv` の残す目標は `e.toEquiv σ` であって `e σ` ではない（2026-09-07、同上）

`e : G ≃* G'` を `Fintype.sum_equiv e.toEquiv f g h` に渡すと、`h` の目標は

```
f σ = g (e.toEquiv σ)
```

になる。ここで `rw [MulEquiv.coe_toEquiv]` を当てると

```
Tactic `rewrite` failed: Did not find an occurrence of the pattern ⇑↑?f
```

（`e.toEquiv σ` は `⇑↑e σ` の形ではなく `Equiv` の直接適用なので当たらない）。
★**直し**は `rw` をやめて `congrArg` で作る:

```lean
(fun σ => (congrArg (fun t => truncENat t (n + 1)) (ramIndex_congr_equiv e ψ hcompat ϖ σ)).symm)
```

★`e.toEquiv σ` と `e σ` は **defeq** なので、補題側の `e σ` の形のまま `congrArg` に渡せば通る。

## #164 `Nat.card (G ⧸ E.fixingSubgroup) = [E : k]` に `restrictNormalHom` は要らない（2026-09-07、同上）

素朴には `AlgEquiv.restrictNormalHom_surjective` + `IntermediateField.restrictNormalHom_ker` +
`QuotientGroup.quotientKerEquivOfSurjective` + `IsGalois.card_aut_eq_finrank` を並べたくなるが、
mathlib には

```lean
IntermediateField.finrank_eq_fixingSubgroup_index   -- FieldTheory/KrullTopology.lean:312
  (L : IntermediateField k K) [IsGalois k K] : Module.finrank k L = L.fixingSubgroup.index
```

が在り、★**`[Normal k L]` を要求しない**（`[IsGalois k K]` だけ）。
`Subgroup.index_eq_card` と繋げば 1 行で `Nat.card (G ⧸ L.fixingSubgroup) = finrank k L` になる。
★★これで **#154（`restrictNormalHom` を含む目標に `exact` を当てると 200000 heartbeats で落ちる）を
丸ごと回避できる** —— `restrictNormalHom` を 1 度も書かない。
★2 層の中間体（`IntermediateField.restrict h`）は `IntermediateField.restrict_algEquiv` +
`LinearEquiv.finrank_eq` で下の層に戻す（#59 の定型 (c)）。

## #165 商群の作用まわりの道具は `Found/PGC/HasseArfInduction.lean` の §4 に**一式ある**（2026-09-07、同上）

★#158（同じ名前で書いて `already been declared` を出させる）で **6 本**が在庫と分かった:

| 在庫 | 内容 |
|---|---|
| `quotientMulSemiringActionOfTrivial` | `H` が `C` に自明に作用するとき `MulSemiringAction (G ⧸ H) C`（`compHom` + `QuotientGroup.lift`） |
| `quotientSMul_mk` / `quotientSMul_mk_fixedRing` | `(σ̄) • c = σ • c`（どちらも `rfl`） |
| `ramIndex_quotient_mk` | `i_ϖ(σ̄) = i_ϖ(σ)` |
| `mem_ramificationGroupReal_quotient_mk_iff` | 下付き分岐群の商版 |
| `herbrandPhiGroup_quotient_eq` | `φ_{G/H} = φ_G`（`ϖ` が固定環の元のとき） |

★★**`MulSemiringAction (G ⧸ H) ↥(fixedRing B H)` は `instance` にしないのが木の流儀**である
（`HasseArfInduction` が「文脈ごとに選ぶものなので `letI` で入れる」と明記している）。
新しい補題も `[MulSemiringAction (G ⧸ H) C]` + `hq : ∀ σ c, (mk σ) • c = σ • c` を
**仮定として受け取る**形に揃えること（`ramIndex_quotient_mk` と同じ形）。

## #166 `refine (MulEquiv.map_eq_one_iff _).mp ?_` は `MulOneClass ?m` で詰まる —— 先に `have` で等式を作る（2026-09-07、`hρ` を外す / Prop 6.14 n=1）

```lean
-- ✗ typeclass instance problem is stuck / MulOneClass ?m.368
refine (MulEquiv.map_eq_one_iff _).mp ?_
rw [← hu]; exact (QuotientGroup.eq_one_iff u).mpr hmM

-- ✓ 等式を先に作ってから当てる
have h1 : ρ σ = 1 := by rw [← hu]; exact (QuotientGroup.eq_one_iff u).mpr hmM
exact (MulEquiv.map_eq_one_iff _).mp h1
```

`MulEquiv.map_eq_one_iff (h : M ≃* N) {x : M} : h x = 1 ↔ x = 1` は、**結論 `x = 1` からは
`N` も `h` も決まらない**。`refine` の `_` はその時点で `?N`・`?h` のままなので
`MulOneClass ?N` の探索が走って止まる。★**逆向き（`.mp` の入力側）を先に `have` で
具体化してから当てる**と一発で通る。同じ形は `MulEquiv.map_eq_one_iff` に限らず、
「片側にしか出てこない型」を持つ iff 全般で起きる。

## #167 `q^a ≤ q^b ↔ a ≤ b` には名前が在る（`Nat.pow_le_pow_iff_right`）—— #161 の警告は `≤` 版には当たらない（2026-09-07、同上）

`#161` は「`a^k = a^n → k = n` は `exact?` が見つけられない」と記録しているが、
★**不等号版はそのまま在る**:

```lean
Nat.pow_le_pow_iff_right : 1 < b → (b ^ m ≤ b ^ n ↔ m ≤ n)
```

★`1 < b` の位置に `hq2 : 2 ≤ pp ^ ff` をそのまま渡せる（ℕ では同じ命題）。
`ℕ∞` に持ち上がっている目標には `Nat.cast_le` を先に当てて ℕ に落とす:

```lean
rw [hram, Nat.cast_le, Nat.pow_le_pow_iff_right hq2]   -- ((q^(i+1):ℕ):ℕ∞) ≤ ((q^i₀:ℕ):ℕ∞) ↦ i+1 ≤ i₀
```

## #168 フィルターの「切れ目」は `∃ i, P i ∧ ¬ P (i+1)` ではなく **`∃ i₀, ∀ j, (P j ↔ j ≤ i₀)`** で出す（2026-09-07、同上）

主単数の列 `u ∈ 1+𝔭^j` のような「1 段ずつ降りる ℕ 上の述語」から段を取り出すとき、
切れ目の存在だけを出すと**消費側で毎回 `j ≤ i₀` を組み直す羽目になる**。
同値の形で出しておくと、そのあとが `exact hiff j` で済む。

```lean
theorem exists_le_iff_of_antitone_natPred (P : ℕ → Prop)
    (hanti : ∀ i : ℕ, P (i + 1) → P i) (h0 : P 0) (m : ℕ) (hm : ¬ P m) :
    ∃ i₀ : ℕ, i₀ < m ∧ ∀ j : ℕ, (P j ↔ j ≤ i₀)
```

★`i₀` の**一意性も同値の形から自動で出る**（2 つあれば互いに `≤`）。
★★実測: `lean_check` **0.07 秒・一発**。`Found/PGC/LubinTateRamificationBookkeeping.lean` の §0。
★ℕ の切り詰め引き算を出さずに `0 ≤ i₀ < m` を扱うには、そのあと
`obtain ⟨k, hk⟩ : ∃ k, M = i₀ + k := ⟨M - i₀, by omega⟩; subst hk` とする
（`M + 1` が**そのまま** `i₀ + k + 1` になるので、段 1 の補題が引数の書き換え無しで当たる）。


## #169 3 層の部分型は「周囲の体への 1 本の `RingHom`」に潰す —— #59 と #69 を**同時に**避ける（2026-09-07、B5 の穴 2 / Thm 6.15）

`↥(fixedRing ↥(adjoinIntegers K x) H)` は `K.closure` から見て**部分型 3 層**である。
このまま `rfl` で層をまたぐと #59（kernel が止まる）、
体の側に出ようとすると #69（`adjoinField`/`adjoinIntegers` の境界）に当たる。

★**逃げ方: 3 つの包含 `RingHom` を合成して 1 本にし、以後は「`K.closure` の元が等しいか」しか問わない。**

```lean
noncomputable def fixedRingToClosure (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier ↥K.carrier⟮x⟯]
    (H : Subgroup (↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)) :
    ↥(fixedRing ↥(adjoinIntegers K x) H) →+* K.closure :=
  (K.carrier⟮x⟯.val.toRingHom).comp
    ((adjoinIntegers K x).subtype.comp (fixedRing ↥(adjoinIntegers K x) H).subtype)
```

★**単射性は `fun _ _ hab => Subtype.ext (Subtype.ext (Subtype.ext hab))` の 1 行**で出る（層の数だけ `Subtype.ext`）。

そのうえで、抽象核

```lean
noncomputable def ringEquivOfRangeEq (f : A →+* L) (g : B →+* L)
    (hf : Function.Injective f) (hg : Function.Injective g)
    (hr : f.range = g.range) : A ≃+* B
```

（`RingEquiv.ofBijective f.rangeRestrict` → `RingEquiv.subringCongr hr` → `(… g …).symm`）に
**像の一致だけ**を渡せば環同型が出る。★可換環しか出てこないので `lean_check` **0.29 秒・一発**。
作用の保存も同じ核で書ける：

```lean
(hact : ∀ (σ : G) (a : A) (b : B), f a = g b → f (σ • a) = g (e σ • b)) →
  ringEquivOfRangeEq f g hf hg hr (σ • a) = e σ • ringEquivOfRangeEq f g hf hg hr a
```

★`A`・`B` の作用は `L` の上の作用ではないので、`hact` の形（「`L` の中で一致する元は、
作用させても `L` の中で一致する」）にするのが要点である。
★実測: `Found/PGC/FixedRingAdjoinIso.lean`。**#59 にも #69 にも 1 度も当たらなかった。**

## #170 `MulEquiv.irreducible_iff (f := (e : _ ≃* _))` の `_ ≃* _` は metavariable になる（2026-09-07、同上）

`RingEquiv` を `MulEquiv` として渡すとき、型を `_` で置くと `f` が決まらず

```
Type mismatch: (MulEquiv.irreducible_iff ?m.110).mpr hα has type Irreducible (?m.110 α₀)
```

になる。★`RingEquiv` は `MulEquivClass` なので**そのまま渡せばよい**：

```lean
(MulEquiv.irreducible_iff (f := (fixedRingAdjoinEquiv K x x₀ h).symm) (x := α₀)).2 hα
```

（在庫 `addVal_ringEquiv` のように `(ψ : C ≃* C')` と**具体的な型**で書くのも可。
`_ ≃* _` だけが駄目である。）

## #171 `AlgEquiv.restrictNormalHom` は `MonoidHom.mk'` に展開されるので `rw` が刺さらない —— `σ.restrictNormal E` で `have` を作る（2026-09-07、同上）

`restrictNormalHom E σ` を含む目標に `AlgEquiv.restrictNormal_commutes` を `rw` しようとすると

```
Did not find an occurrence of the pattern ↑((σ.restrictNormal ↥(restrict h)) …)
in the target expression … ((MonoidHom.mk' (fun χ => χ.restrictNormal ↥(restrict h)) ⋯) σ) …
```

で落ちる。★**`restrictNormalHom E σ = σ.restrictNormal E` は `rfl`** なので、
最初から `σ.restrictNormal E` の形で `have … := rfl` を作り、それを `rw` すればよい。

```lean
have hstep : galQuotientEquiv h (QuotientGroup.mk σ) y
    = (IntermediateField.restrict_algEquiv h).symm
      (σ.restrictNormal ↥(IntermediateField.restrict h)
        (IntermediateField.restrict_algEquiv h y)) := rfl
```

★★**#154（`restrictNormalHom` に `exact` を当てると 200000 heartbeats で落ちる）は
`K₁ := K.closure`（無限次元）のときの話である。** `K₁ := ↥K.carrier⟮x⟯`（**有限次拡大**）では
起きなかった（`galQuotientEquiv` の定義と `coe_galQuotientEquiv_apply` で 0.74 秒）。
★#154 を見て `restrictNormalHom` を避ける前に、**大きい方の体が有限次かどうかを見ること。**

## #172 像の等式で作った環同型に `exact` するときは `f`・`g` を**明示する**（2026-09-07、同上）

`apply_ringEquivOfRangeEq _ _ _ _ _ c` を、結論を**展開した coercion の形**で書いた補題に
当てると metavariable が決まらない：

```
apply_ringEquivOfRangeEq ?m.115 ?m.116 … has type ?m.116 (… c) = ?m.115 c
but is expected to have type ↑↑((fixedRingAdjoinEquiv K x x₀ h) c) = ↑↑↑c
```

★`f`・`g`（と単射性・像の等式）を**全部書く**と一発で通る。定義的には等しいので `rfl` 系の
苦労は要らない——**足りないのは型推論の手がかりだけ**である。


## #173 「第 1 引数が explicit」で顔がまったく違うエラーになる 3 つ（2026-09-07、B6 / Thm 6.15）

`Subgroup.Normal.of_commutator_le` は **群 `G` が第 1 explicit 引数**である。
`Subgroup.Normal.of_commutator_le h` と書くと

```
Application type mismatch: The argument h has type
  commutator Gal(L/k) ≤ E.fixingSubgroup of sort `Prop`
  but is expected to have type Type ?u.29
```

★「Prop を Type の位置に置いた」という顔になるので、**`h` の型が悪いのだと誤読しやすい**。
正しくは `Subgroup.Normal.of_commutator_le _ h`。

`IntermediateField.lift_adjoin_simple` / `IntermediateField.lift_top` も
**底体 `F` が第 1 explicit 引数**である。`lift_adjoin_simple M α` と書くと

```
the argument α has type ↥M but is expected to have type IntermediateField ↥M ?m.44
```

★`M` が `F` の位置に吸われて `E := ↥M` と推論されるためで、
「α が中間体でない」という**無関係な顔**になる。正しくは
`IntermediateField.lift_adjoin_simple k M α` / `IntermediateField.lift_top k M`。

★★**共通の直し方は `#check @名前` を 1 回叩くこと**（0.01 秒）。
エラー文から explicit/implicit を推測しない。★3 件とも 1 往復で片付いた。

## #174 「有限次拡大は単項生成」は `lift` に閉じ込めれば 3 行（#59 の定型 (e)）（2026-09-07、同上）

mathlib の `Field.exists_primitive_element k ↥M` は `∃ α : ↥M, k⟮α⟯ = ⊤` という
**`M` の中の**主張なので、そのままでは周囲の体 `L` の中の
`IntermediateField.adjoin k {y} = M` にならない。`lift` で降ろす:

```lean
theorem exists_adjoin_singleton_eq {k L : Type*} [Field k] [Field L] [Algebra k L]
    (M : IntermediateField k L) [FiniteDimensional k ↥M] [Algebra.IsSeparable k ↥M] :
    ∃ y : L, IntermediateField.adjoin k ({y} : Set L) = M := by
  obtain ⟨α, hα⟩ := Field.exists_primitive_element k ↥M
  exact ⟨(α : L), by
    rw [← IntermediateField.lift_adjoin_simple k M α, hα, IntermediateField.lift_top k M]⟩
```

★2 層をまたぐ `rfl` を 1 つも書かないので **#59 に触れない**（0.15 秒）。
★`Algebra.IsSeparable k ↥M` は `[IsGalois k L]` があれば `inferInstance` で出る。

## #175 下付きの「for a large m」から上付きへは 6 行で乗り換えられる（2026-09-07、同上）

原典が「`Gal(K′/L)_m = {id}` for a large m」と書くとき、消費側（Cor 6.13(ii)）が
要求するのは**上付き** `G^m` である。木にあるのは下付きの
`exists_lowerRamificationGroup_eq_bot` だけだが、**新しい数学は要らない**:

1. `G_N = ⊥` なる `N` を取る。
2. `upperRamificationGroup_herbrandPhiGroup` で `G^{φ_G(N)} = ramificationGroupReal α N`。
3. `ramificationGroupReal_eq_of_mem_Ioc` で実数添字を `ℕ` に戻して `= G_N = ⊥`。
4. `exists_nat_ge` で `φ_G(N) ≤ k+1` なる `k` を取り、`upperRamificationGroup_antitone`。

★`lean_check` 0.23 秒・一発。★結論を `((k+1 : ℕ) : ℝ)` の形にしておくと
B5 系の消費側にそのまま入る（`ℕ` の切り詰め引き算を書かない、#102）。

## #176 「完備化から代数へ降ろす」に Krasner／Lemma 2.2(iii) は要らない —— 無限次 Galois で 5 行（2026-09-07、Y23 / Yoshida Def 4.10）

原典が「`Ê ∩ K^sep = E` だから完備側の一意性が代数側に降りる」と書く箇所は、
henselian 体上の Krasner の補題（重い。mathlib の `IsKrasner` は `[CompleteSpace K]` を
要求するので、完備でない `K^ur` を底にするとそのままでは使えない）を経由するように見える。
★**経由しなくてよい。**

`E, E′ ⊆ Ω`（`Ω/k` は無限次 Galois）で「像の生成する完備体が一致する」ことから
`E′ ⊆ E` を出すには、**`E` を固定する `σ ∈ Gal(Ω/k)` を完備化へ延長する**だけでよい:

1. `σ` はスペクトルノルムで等長（`spectralNorm_eq_of_equiv`）→
   `UniformSpace.Completion.mapRingEquiv` で完備化へ延長。
2. `σ` が稠密な部分体を各点固定 → 延長はその完備化を各点固定（`DenseRange.equalizer`）。
   ここで `AlgEquiv.ofRingEquiv` を噛ませると**完備な底体上の `AlgEquiv`** になる。
3. 生成元も固定なら `IntermediateField.adjoin_le_iff` で `adjoin` 全体を固定（下の抽象核）。
4. 完備化への埋め込みは単射なので `σ x = x` に戻る。
5. `InfiniteGalois.fixedField_fixingSubgroup` で `x ∈ E`。

★抽象核は 2 本だけで、どちらも分岐・付値の語彙が 0 語:

```lean
-- 1 つの AlgEquiv の固定体（mathlib の fixedField は Subgroup を取るので 1 元版が要る）
def fixedIntermediateField (φ : Ω ≃ₐ[k] Ω) : IntermediateField k Ω  -- carrier := {x | φ x = x}
theorem algEquiv_eq_self_of_mem_adjoin (φ) (S) (hS : ∀ s ∈ S, φ s = s) :
    x ∈ IntermediateField.adjoin k S → φ x = x
theorem mem_of_forall_fixingSubgroup_eq [IsGalois k Ω] (E) (x)
    (h : ∀ σ : Ω ≃ₐ[k] Ω, (∀ y ∈ E, σ y = y) → σ x = x) : x ∈ E
```

★実測（Y23）: 抽象核 2 本で **0.13 秒・ほぼ一発**、降下の具体層 **0.65 秒・一発**、
ファイル全体 529 行が `leanfile.mjs` で**一発 ok**。原典が Krasner で畳んだ段が
**新しい数学ゼロ**で済んだ。

★付随する小さな穴 2 つ:

* `SetLike` の `≤` は**関数ではない**ので `le_sup_right hy` は `Function expected` になる。
  `SetLike.le_def.mp le_sup_right hy` と書く（2342 行の「型注釈が要る」と同じ穴の別の顔）。
* `IntermediateField` の構造体フィールドに `neg_mem'` は**無い**（`inv_mem'` はある）。
  `where` で書くと `neg_mem' is not a field of structure IntermediateField` が出る。

## #177 `rw [← 補題] at h` が「型クラス合成に失敗」の顔で落ちる —— 逆向きは**明示引数を全部渡す**（2026-09-07、Λ8 / Yoshida Def 4.10）

**症状**: 次の 2 行が並んでいる。

```lean
rw [range_abelianGalEquivProd_artinMap] at hmem      -- 正順: 通る
rw [← range_abelianGalEquivProd_artinMap] at hmem    -- 逆順: 落ちる
```

出るエラーは `rewrite failed` ではなく

```
error: failed to synthesize
  ExpChar (IsLocalRing.ResidueField ↑𝒪[K.carrier]) ?pp
```

**原因**: 逆向きの `rw` は補題の**右辺**とパターンマッチする。右辺が
`Set.univ ×ˢ Set.range (fun n : ℤ => arithFrobenius K ^ n)` のように
**暗黙引数（`pp` / `ff` / `hq` …）を 1 つも含まない**と、それらが metavariable
のまま残り、`[ExpChar … pp]` の合成が `?pp` に対して走って落ちる。
★**「rw が失敗した」という顔をしない**のが厄介で、インスタンスの不足を疑って
`haveI` を足しても直らない。

**直し方**: 逆向きのときは補題に**明示引数を全部渡す**。

```lean
rw [← range_abelianGalEquivProd_artinMap K hq hπmax hπne0 f hf0 hf1 hf hM] at hmem
```

★一般則: **右辺に現れない暗黙引数を持つ補題は、`←` で使うときだけ明示化が要る。**
正順で通ったからといって逆順が通るとは限らない（左辺には出ていた、というだけ）。
★同型の `.injective.eq_iff` を `←` で使うとき（`rw [← e.injective.eq_iff]`）も同じ穴に
落ちうるが、そちらは `e` を明示的に書いているので実際には落ちなかった。

## #178 `ring` は `CommMonoid` の乗法だけの等式では `ring_nf made no progress` になる（2026-09-07、Y24 / Yoshida Cor 4.9 後半）

**症状**: 抽象核を `[CommMonoid S]` で書くと、`a * (b * c) = b * (a * c)` のような
**乗法だけ**の等式に `ring` を当てても

```
error: `ring_nf` made no progress on the goal
```

が出る。★「`ring` が失敗した」ではなく「進まなかった」という顔なので、
仮定が足りないのか式が違うのかの区別がつかない。

**原因**: `ring` は `CommSemiring`（少なくとも加法つき）を要求する。
`CommMonoid` には `+` が無いので `ring` の正規化器が起動しない。

**直し方**: 名前で書く。3 つで足りる。

```lean
calc c * (x * t) = x * (c * t) := mul_left_comm _ _ _
  _ = x * (t' * θ)             := by rw [hθ]
  _ = x * t' * θ               := (mul_assoc _ _ _).symm
```

★`mul_comm` / `mul_left_comm` / `mul_assoc` で、実測では 3 行で閉じた。
★関連: `mul_right_cancel h` は `h : a * c = b * c` の形しか受けない。
`h : a * c = c * b` のときは `mul_right_cancel (h.trans (mul_comm _ _))` と書く
（これも「Application type mismatch」という無関係な顔のエラーになる）。

## #179 抽象核を `MulAction S N` で書くと消費側（冪級数）がインスタンスを持たない —— 合成則版を**併置**する（2026-09-07、同上）

**症状**: 「`[a] ∘ [b] = [ab]`」を `MulAction S N` として抽象化すると綺麗に書けるが、
消費側の `[a]` は `PowerSeries` で、`MulAction` インスタンスは木にもmathlib にも無い。
`Function.End N` を経由しようとしても **`Function.End N` は `CommMonoid` ではない**ので
抽象核の `[CommMonoid S]` に合わない（`MulAction.compHom` は使えるが、
`letI` は宣言の境界を越えない（#133）ので statement には書けない）。

**直し方**: 抽象核を**2 本**置く。中身は同じ 1 行の計算である。

```lean
-- (i) MulAction 版（他の抽象核と繋ぐとき用）
theorem smul_spec_transfer {S} [CommMonoid S] {N} [MulAction S N] (σ : N → N) … 

-- (ii) 合成則版（消費側が実際に使う）
theorem map_spec_transfer {S} [CommMonoid S] {N}
    (br : S → N → N) (hbr : ∀ a b n, br a (br b n) = br (a * b) n) (σ : N → N) … 
```

★(ii) は `br 1 = id` も `br` の単射性も使わない。`MulAction` の公理のうち
**合成則しか要らない**ことが証明を書いてみると分かる。
★実測（Y24）: (i) と (ii) の証明はどちらも `rintro` + `rw` の 2 行、
`lean_check` 0.18 秒で一発。**併置のコストは 15 行**。


## #180 section の `variable` 仮説は `def` では拾われ `theorem` では拾われない —— 顔は `Unknown identifier`（2026-09-07、Λ9 / pGC Prop 1.1）

**症状**: 同じ section の中で

```lean
variable {G : Type*} [CommGroup G] (e : G ≃* (U × Z)) {m : ℕ} (hm : m ≠ 0)

noncomputable def torsionEquivUnitsTorsion : … :=          -- ★通る
  (powTorsionCongr e m).trans (powTorsionProdEquiv m (fun _ hb => zhat_eq_one … hm hb))

theorem snd_eq_one_of_mem_torsion (x : ↥(powTorsion G m)) : (e (x : G)).2 = 1 :=
  zhat_eq_one_of_pow_eq_one hm …                            -- ★`Unknown identifier hm`
```

★**`def` は本体を走査して section 変数を含めるが、`theorem` は結論（型）しか見ない。**
`hm` は結論に現れないので落ちる。#151（「section の `variable` 仮説は黙って落ち、
別の顔で出る」）の**具体的な境界**がこれである。

**直し方**: `include hm in` を置く。★**`include … in` は docstring の前**
（#107/#154 と同じ位置。docstring の後ろに書くと
`unexpected token 'include'; expected 'lemma'` になる）。

```lean
include hm in
/-- 捩れ元の `Ẑ` 成分は自明。 -/
theorem snd_eq_one_of_mem_torsion … 
```

★実測（Λ9）: 531 行のファイルで踏んだのはこの 1 箇所だけだった。
`def` で通っていた直後に同じ変数を使う `theorem` を書くと必ず踏む。

## #181 部分型の `.2` は `x ∈ ker f` であって `↑x ^ m = 1` ではない —— `rw` は defeq を見ない（2026-09-07、同上）

**症状**: `powTorsion A m := (powMonoidHom m : A →* A).ker` と定義したとき、
`x : ↥(powTorsion A m)` の `x.2` は **`(powMonoidHom m) ↑x = 1`** という形をしている。
`↑x ^ m = 1` とは defeq だが**構文的には別物**なので

```lean
have h : e (x : G) ^ m = 1 := by rw [← map_pow, x.2, map_one]
```

が

```
Tactic `rewrite` failed: Did not find an occurrence of the pattern
  (powMonoidHom m) ↑x
in the target expression
  e (↑x ^ m) = 1
```

で落ちる（★`rw` は `x.2` の**左辺**を探しに行く）。

**直し方**: 一度 `have` で受け直すと defeq が効く。

```lean
have hx : (x : G) ^ m = 1 := x.2      -- ★これは通る（defeq でチェックされる）
have h : e (x : G) ^ m = 1 := by rw [← map_pow, hx, map_one]
```

★`exact`・`apply`・関数適用（`hΨ _ x.2`）では `x.2` をそのまま渡してよい。
落ちるのは `rw` / `simp only` のように**構文照合する**タクティクだけである。
★#160（`(ϖ.2 : ‖ϖ‖ ≤ 1)` は不可）と裏表の関係にある。

## #182 `where` の構造体インスタンス記法は、先に書いたフィールドが後続の期待型に伝わらないことがある —— `by exact` で包む（2026-09-07、同上）

**症状**: `ContinuousMulEquiv.prodCongr` を作ろうとして

```lean
def prodCongr (e : A ≃ₜ* C) (f : B ≃ₜ* D) : (A × B) ≃ₜ* (C × D) where
  toMulEquiv := e.toMulEquiv.prodCongr f.toMulEquiv
  continuous_toFun := (e.continuous_toFun.comp continuous_fst).prodMk … 
```

と書くと

```
Type mismatch … but is expected to have type
  Continuous (MulEquiv.prodCongr ?m.28 ?m.35).toFun
```

★**`?m` が残っている**。`toMulEquiv` を先に書いたのに、後続フィールドの期待型は
まだ metavariable のままで、`MulEquiv.prodCongr` の暗黙引数が決まっていない。
`{ … with … }` 記法に替えても同じで、さらに `MulOneClass A` の合成失敗まで増える。

**直し方**: 後続フィールドを **`by exact …`** で包む（タクティクブロックは
構造体の全フィールドが決まったあとに走るので `?m` が解ける）。

```lean
def prodCongr (e : A ≃ₜ* C) (f : B ≃ₜ* D) : (A × B) ≃ₜ* (C × D) where
  toMulEquiv := e.toMulEquiv.prodCongr f.toMulEquiv
  continuous_toFun := by
    exact (e.continuous_toFun.comp continuous_fst).prodMk (f.continuous_toFun.comp continuous_snd)
  continuous_invFun := by
    exact (e.symm.continuous_toFun.comp continuous_fst).prodMk
      (f.symm.continuous_toFun.comp continuous_snd)
```

★ついでの在庫測定: **`MulEquiv.prodCongr` は `MulOneClass` を要求する**（`Mul` では足りない）。
`ContinuousMulEquiv.prodCongr` は **mathlib に無い**（`#158` 空振り、
`grep -n "ContinuousMulEquiv" .cache/mathlib-index.txt` に `prodCongr` の行が無い）。

## #183 `Gal(K^ab/K)` の `CommGroup` は `open scoped IsMulCommutative` で出る —— 自作しない（2026-09-07、同上）

**症状になりかけたこと**: `powTorsion A m` は `[CommGroup A]` を要求するが、
`Gal(K^ab/K) = (abelianClosure K ≃ₐ[K.carrier] abelianClosure K)` に付いているのは
`Group` と、木が用意した **`IsAbelianGalois K.carrier (abelianClosure K)`** だけである。
「`CommGroup` インスタンスが無いから `mul_comm` を仮説で受け取ろう」と考えたくなる。

**実測**: mathlib の `Algebra/Group/Defs.lean:1409` に

```lean
scoped instance (priority := 50) {G : Type*} [Group G] [IsMulCommutative G] : CommGroup G
```

が在る（`IsAbelianGalois` は `IsMulCommutative Gal(L/K)` を extends している）。
★**`open scoped IsMulCommutative` を 1 行足すだけで `CommGroup Gal(K^ab/K)` が出る。**

引き当ては名前空間 1 回 grep（#117(ii)）:

```
grep -n "IsMulCommutative" .cache/mathlib-index.txt | head
```

★これで `IsAbelianGalois` の定義行（`extends IsGalois K L, IsMulCommutative Gal(L/K)`）と
scoped instance の行が同時に出る。**「無い」と書く前に必ずこの 1 手を踏むこと。**

## #184 `PowerSeries` のままだと代入の結合律が `∀ x y z` の形で立たない —— 定数項 `0` の部分型へ移すと純抽象核が当たる（2026-09-07、Y26 = Cor 4.9 具体化 1 本目）

**症状**: 「合成の結合律だけで出る」主張（絡み作用素の合成・左簡約など）を
`comp : M → M → M` + `hassoc : ∀ x y z, comp (comp x y) z = comp x (comp y z)`
という**純抽象核**（型 1 つと二項演算 1 つ、`#print axioms` が
`does not depend on any axioms`）に切り出したい。ところが
`M := PowerSeries A`・`comp x y := PowerSeries.subst y x` と置くと `hassoc` が**立たない**:
mathlib の

```lean
PowerSeries.subst_comp_subst_apply (ha : HasSubst a) (hb : HasSubst b) (f) :
    subst b (subst a f) = subst (subst b a) f
```

は `HasSubst`（= 定数項が冪零）を**両方に**要求するので `∀ x y z` にならない。

**直し方**: 部分型へ移す。`HasSubst` が自動になり `hassoc` が無条件になる。

```lean
def SubstNilp (A : Type*) [CommRing A] := {p : PowerSeries A // PowerSeries.constantCoeff p = 0}

noncomputable def substComp (x y : SubstNilp A) : SubstNilp A :=
  ⟨PowerSeries.subst y.1 x.1, PowerSeries.constantCoeff_subst_eq_zero y.2 x.1 x.2⟩

theorem substComp_assoc (x y z : SubstNilp A) :
    substComp (substComp x y) z = substComp x (substComp y z) :=
  Subtype.ext (PowerSeries.subst_comp_subst_apply (hasSubst_of_constantCoeff_zero y.2)
    (hasSubst_of_constantCoeff_zero z.2) x.1)

noncomputable def substOne : SubstNilp A := ⟨PowerSeries.X, PowerSeries.constantCoeff_X⟩
theorem substOne_comp (x : SubstNilp A) : substComp substOne x = x :=
  Subtype.ext (PowerSeries.subst_X (hasSubst_of_constantCoeff_zero x.2))
```

★戻すのは `congrArg Subtype.val (核 … ⟨f, hf0⟩ ⟨g, hg0⟩ … (Subtype.ext h₁) …)` の 1 行。
`(substComp x y).1 = subst y.1 x.1` は `rfl` なので、入れ子も含めて何も展開しなくてよい。
★**部分型は 1 層なので #59（中間体 2 層の `rfl`）には当たらない**。実測 0.22 秒で一発。
★副産物: 結合律を使う場所がファイル中 `substComp_assoc` の 1 箇所に閉じる。
以後の証明に `subst_comp_subst_apply` は 1 度も現れない。
★**`noncomputable` を忘れると `depends on 'PowerSeries.X', which is 'noncomputable'`**
という顔で `substOne` の側が落ちる（結合律とは無関係な行が赤くなる）。

## #185 Git Bash の `grep` / `cat -A` は CR を落とす —— CRLF の判定に使ってはならない（2026-09-07、同上）

**症状**: `lean/ABC3/Found.lean` は CRLF（#141）。import を 1 行足す前に
行末を確かめようとして

```
$ sed -n '1799,1800p' ABC3/Found.lean | cat -A
import ABC3.Found.PGC.LubinTateReciprocityIndependence$
import ABC3.Found.PGC.CyclotomicFromAbelianization$
```

と出た（`^M$` ではない）。**「最近足された import だけ LF なのだ」と誤読して
`
` で挿入しようとし、Python の `assert b.count(anchor)==1` が `0` で落ちた。**

**原因**: この環境の Git Bash の `grep` / `cat -A` はパイプの途中で CR を捨てる。
`sed -n 'Np' | cat -A` も同じ。★行末は**バイトで見ないと分からない**。

**直し方**: Python でバイト列を直接数える。

```python
b = open('ABC3/Found.lean','rb').read()
print('crlf', b.count(b'

'), 'lf', b.count(b'
'))   # 1806 1806 ⇒ 全部 CRLF
```

★挿入も同じスクリプトの中で `anchor = b'import …

'` として行う
（`sed -i` は Windows 側で行末を書き換えることがある）。
★挿入後に `crlf == lf` が保たれていることを**必ず測り直す**（実測: 1806/1806 → 1807/1807）。

## #186 `Subsingleton (ZMod (p ^ 0))` はインスタンス探索が出さない —— `rw [pow_zero]` してから `infer_instance`（2026-09-07、Λ11 / pGC Prop 1.1）

**症状**: `n = 0` の場合分けで `ZMod (p ^ 0)` の元の等式を潰そうとして

```lean
rcases Nat.eq_zero_or_pos n with rfl | hn
· exact Subsingleton.elim _ _
```

と書いたら `failed to synthesize Subsingleton (ZMod (p ^ 0))`。

**原因**: インスタンス探索は `p ^ 0` を `1` に**簡約しない**（`Nat.pow` は
インスタンス合成の中では簡約されない）。`ZMod 1` の `Subsingleton` は在るのに引けない。

**直し方**: 先に書き換えてから合成する。

```lean
· haveI : Subsingleton (ZMod (p ^ 0)) := by rw [pow_zero]; infer_instance
  exact Subsingleton.elim _ _
```

★同型の顔は `Fintype (ZMod (p ^ 0))` や `NeZero (p ^ 0)` でも出る。
`NeZero` は `⟨pow_ne_zero _ hp.ne_zero⟩` で作る（こちらは `pow_ne_zero` が在る）。

## #187 mathlib の**インスタンス**の import 漏れは「型クラス合成に失敗」の顔で出る（#68 のインスタンス版）（2026-09-07、同上）

**症状**: `InfiniteGalois.normalAutEquivQuotient`（閉正規部分群による商 ≅ 固定体の Galois 群）に
交換子群の位相的閉包を渡したら

```
failed to synthesize instance of type class
  (↑{ toSubgroup := (commutator Gal(E/k)).topologicalClosure, isClosed' := ⋯ }).Normal
```

`Subgroup.Normal` の合成が落ちているように見えるので「正規性を自分で証明しなければ」と
読みたくなる。**そうではない。**

**原因**: `instNormalCommutatorClosure : (commutator G).topologicalClosure.Normal` は
`Mathlib/Topology/Algebra/Group/TopologicalAbelianization.lean` にあり、
**その 1 ファイルを import していないだけ**だった（#68 の「Unknown constant は import 漏れ」の
**インスタンス版**。定数名が出ないぶん気づきにくい）。

**直し方**: 名前空間 1 回 grep で担い手を探す。

```
$ grep -n "instNormalCommutatorClosure\|TopologicalAbelianization" .cache/mathlib-index.txt
3194:abbrev   TopologicalAbelianization   Topology/Algebra/Group/TopologicalAbelianization.lean:42
61806:instance instNormalCommutatorClosure  .../TopologicalAbelianization.lean:37
```

★共有 REPL は他の agent の import を持っているだけなので、
**「REPL で合成できない」は「mathlib に無い」ではない**。
REPL では `[N.Normal]` を仮引数に足して抽象核だけ通し、
具体層は `leanfile.mjs`（正しい import つき）で測る——本ノードはこれで 0.06 秒 + 11.4 秒で済んだ。

## #188 `rw … at h` が「motive が型付かない」で落ちる —— 部分型は先に `rintro ⟨a, ha⟩` で潰す（2026-09-07、Λ11 / pGC Prop 1.1）

**症状**: `x : ↥(cyclotome ↥S (p ^ 0))` について

```lean
have ha : (a : TopologicalAbelianization ↥S) ^ (p ^ 0) = 1 := a.2
rw [pow_zero, pow_one] at ha        -- motive is not type correct
```

`↑a` has type `↥(cyclotome (↥S) (p ^ 0))` but is expected to have type
`↥(cyclotome (↥S) _a)` と言われる。

**原因**: `rw` は `p ^ 0` を**すべての出現箇所で**抽象化する。`↑a` の
コアーション関数（`Subtype.val`）の型に `p ^ 0` が入っているので、抽象化すると
`↑a` が型付かなくなる。★#181 の親戚だが顔が違う（あちらは defeq、こちらは motive）。

**直し方**: 部分型を**先に壊す**。台の元の型には `p ^ 0` が出てこない。

```lean
rintro ⟨a, ha⟩ ⟨b, hb⟩
have ha' : a ^ (p ^ 0) = 1 := ha    -- a : TopologicalAbelianization ↥S（p^0 を含まない）
rw [pow_zero, pow_one] at ha'       -- 通る
exact Subtype.ext (ha'.trans hb'.symm)
```

★同じ理由で `rw [Subtype.ext h]`（`x = 1` で `x` を書き換える）も落ちる。
`Subsingleton` インスタンスを作って `Subsingleton.elim` に逃がすのが速い。

## #189 `linear_combination` の残差に「× 2」が出たら、その係数の符号が逆（2026-09-07、Yoshida08 Prop 3.5 ねじれ版）

**症状**: `linear_combination hobs' + π' * hc + ϕ c * hπpow` が

```
ring failed, ring expressions not equal
⊢ π' * π' ^ n * ϕ c * s * s ^ n * 2 - ϕ c * π * π ^ n * 2 = 0
```

**読み方**: 残差の各項に **ちょうど 2 が掛かっている**ときは、
「係数が 0 であるべき」ではなく「**符号が逆**」である（正しい寄与 `-x` を `+x` と
書いたので `2x` が残る）。`ϕ c * hπpow` を `- ϕ c * hπpow` に変えて一発で通った。
★残差を睨んで項を足すより、**2 の有無を見る**ほうが速い。

**あわせて**: `ring` は加群（`Module R M`）の目標では `ring_nf made no progress` で
落ちる（環ではないので当然だが顔が分かりにくい）。`x = (id - T) x + T x` の類は
`abel` を使う。

## #190 `Ideal.map_le_iff_le_comap.mpr (fun y hy => …)` は項で書くと ?m が決まらない（2026-09-07、同）

**症状**:

```lean
exact Ideal.pow_right_mono (Ideal.map_le_iff_le_comap.mpr fun y hy => hϕloc y hy) k
-- Type mismatch: hϕloc y hy : ϕ y ∈ maximalIdeal A
--   but is expected to have type  y ∈ Ideal.comap ?m.241 ?m.243
```

**原因**: `comap f I` の `f` と `I` が決まらないうちに本体を型検査するので、
`mem_comap` の defeq が使えない。

**直し方**: いったん `have` に切り出して **`rw` で `comap` の形にしてから `intro`**:

```lean
have hmapm : Ideal.map ϕ (IsLocalRing.maximalIdeal A) ≤ IsLocalRing.maximalIdeal A := by
  rw [Ideal.map_le_iff_le_comap]
  intro y hy
  exact hϕloc y hy          -- ここでは comap の membership が defeq で通る
```

★ついでに: イデアルの冪の単調性は `Ideal.pow_mono` **ではなく**
`Ideal.pow_right_mono`（`exact?` が 2.2 秒で出した）。

## #191 `show ... from rfl` を三重強制の途中に置くと「`Monoid ?m` で instance が stuck」（2026-09-07、Λ12 = Artin 写像の同変性）

**症状**: `rootsOfUnity` のような **`Subgroup Rˣ`**（= 部分型の中の単数の中の環の元、
強制が 3 段）を扱うとき、次のように書くと落ちる。

```lean
rw [show ((conj x : ↥(rootsOfUnity m ↥L)) : (↥L)ˣ)
      = (MulEquiv.restrictRootsOfUnity σ m x : ↥(rootsOfUnity m ↥L)) from rfl]
-- typeclass instance problem is stuck
--   Monoid ?m.92
```

★エラーは `rfl` が偽であることを言っていない。**`show` の型を elaborate する途中で
強制の 1 段目（`Monoid` の何か）が metavariable のまま**になるだけである。

**直し方**: **強制を最後まで書き下した `@[simp]` 補題を 1 本立てて、`rw` で使う。**

```lean
@[simp] theorem conj_coe (…) :
    (((conj x : ↥(rootsOfUnity m ↥L)) : (↥L)ˣ) : ↥L)      -- ★3 段とも書く
      = σ (((x : (↥L)ˣ) : ↥L)) :=
  MulEquiv.restrictRootsOfUnity_coe_apply _ _            -- mathlib がちょうど持っている
```

そのうえで本体は `refine Subtype.ext (Units.ext ?_)` で 2 段はがしてから
`rw [conj_coe, …]`。★実測: 3 箇所すべてこの形で消えた（往復 1 回）。

★**教訓**: `Subgroup Mˣ` を相手にするときは、
**最初に「一番下まで降ろした coe 補題」を作る**。それ以降の証明は `rw` の直線になる。

## #192 `MulEquiv.trans` は `rw` で展開されない —— 残る goal は `f (e x) = f ((e.trans f') x)`（2026-09-07、同上）

**症状**: `e : A ≃* B`、`f : B ≃* C` について
`(e.trans f) x` を含む goal を `rw` で潰したあと、

```
⊢ conj (e' (E x)) = conj ((E.trans e') x)
```

が残る（`unsolved goals`）。★`rw` は `MulEquiv.trans_apply` を自動では使わない。

**直し方（どちらでもよい）**:

1. 末尾に **`rfl` を 1 行足す**（`trans` の適用は defeq）。
2. ★**`show` で両辺を `trans` を使わずに書き直してから `rw`**。
   ゴールが `∃ e, ∀ g x, …` の形で `e := E.symm.trans e'` のように
   **自分で作った合成を入れている**ときは、こちらでないと `rw` の
   パターンが当たらない（`rewrite failed: did not find an occurrence`）。

```lean
show e' (E.symm (F g x)) = conj (e' (E.symm x))   -- ★trans を消して書く
rw [hx, he g]
```

★実測: 「`∃ e, ∀ …` の `e` に合成を入れる」形の証明を 4 本書いて、
**3 本で (2) が必要**だった。(1) で済むのは goal の左辺だけに `trans` が出る場合。

## #191 抽象核の「イデアル引数」を `_` のまま呼ぶと `whnf` が heartbeat を焼く（2026-09-07、Yoshida08 Lemma 3.4 完備版）

**症状**: 抽象核を具体層から呼ぶだけの 1 行が
`(deterministic) timeout at 'whnf', maximum number of heartbeats (200000)` で落ちる。

```lean
-- 落ちる（`_` は I = IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K)）
exact semilinear_solution_unique_of_isHausdorff _ (unramGalCompletionInt K σ : _ →+* _)
  (fun x hx => maximalIdeal_map_mem_of_ringEquiv (unramGalCompletionInt K σ) x hx) hu h
```

**原因**: `_` に入る `I` は `[IsHausdorff I A]` のインスタンス探索の**索引**でもある。
`A` が `↥(何かの def)` のとき、`I` が未定のままインスタンス探索が始まり、
`unramifiedCompletionInt` の `def` を `whnf` で剥がしにいって帰ってこない。
★同じファイルの `hsolve` 側（`I` を露出しない補題を経由した）は **13 秒で通っていた**ので、
遅いのは環そのものではなく「`_` にした穴」である。

**直し方**: ★**イデアルを露出しない「環同型版」の抽象補題を 1 本挟む**。
穴を `_` にせず、補題の側で `IsLocalRing.maximalIdeal A` に**固定**してしまう。

```lean
theorem semilinear_solution_unique_of_isHausdorff_ringEquiv {A : Type*} [CommRing A]
    [IsLocalRing A] [IsHausdorff (IsLocalRing.maximalIdeal A) A] (e : A ≃+* A)
    {u : A} (hu : u ∈ IsLocalRing.maximalIdeal A) {c c' : A}
    (h : c - u * e c = c' - u * e c') : c = c' :=
  semilinear_solution_unique_of_isHausdorff (IsLocalRing.maximalIdeal A) (e : A →+* A)
    (fun x hx => maximalIdeal_map_mem_of_ringEquiv e x hx) hu h
-- 具体層は `exact … _ringEquiv (unramGalCompletionInt K σ) hu h` の 1 行になり、即座に通る
```

★#150（抽象核の穴を暗黙引数にすると別の顔のエラーが 5 個出る）と同じ教訓の
**インスタンス索引版**である。「明示引数にする」だけでなく
**「インスタンスの索引になる引数は呼び出し側で `_` にしない」**まで言える。

## #192 局所環の `ϕ(𝔪) ⊆ 𝔪` は付値を経由せず `isLocalHom_equiv` で出る（2026-09-07、同）

半線型方程式・Frobenius ねじれで必要になる `ϕ(𝔪) ⊆ 𝔪` は、
`ϕ` が**環同型**でありさえすれば mathlib の instance
`isLocalHom_equiv`（`Mathlib/Algebra/Group/Units/Equiv.lean:214`）だけで出る。
付値・不分岐性・Galois の語彙は 1 つも要らない。

```lean
theorem maximalIdeal_map_mem_of_ringEquiv {A : Type*} [CommRing A] [IsLocalRing A]
    (e : A ≃+* A) (x : A) (hx : x ∈ IsLocalRing.maximalIdeal A) :
    e x ∈ IsLocalRing.maximalIdeal A := by
  rw [IsLocalRing.mem_maximalIdeal, mem_nonunits_iff] at hx ⊢
  exact fun hu => hx (IsLocalHom.map_nonunit (f := e) x hu)
```

★`#print axioms` が `[propext, Quot.sound]`（`Classical.choice` すら要らない）。
★Galois 群の元に対して `σπ = π` の類の性質を探しに行く前に、まずこれを試すこと。


## #193 `(ZMod m)ˣ` に位相インスタンスは無い —— 「連続指標」は**核が開**で書く（2026-09-07、Λ13 = pGC Prop 1.1 の局所 Tate 双対性の道）

副有限群 `Γ` から有限離散群への「連続指標」を Lean で書こうとして

```lean
example (K : PAdicLocalField p) (n : ℕ) : Prop := Continuous (cycloCharUnitsModPow K n)
-- failed to synthesize instance of type class
--   TopologicalSpace (ZMod (p ^ n))ˣ
```

で止まる。★`ZMod m` にも `(ZMod m)ˣ` にも位相インスタンスは**無い**（2026-09-07 実測）。
自分で `⊥`（離散）インスタンスを入れると木全体に漏れるので、やってはならない。

★**逃げ道**: 離散有限群への準同型については「連続」と「核が開」は同値なので、
**`IsOpen (ψ.ker : Set G)` と書く**。位相インスタンスは `G` 側にしか要らない。

```lean
theorem isOpen_ker_cycloCharHom (K : PAdicLocalField p) (n : ℕ) :
    IsOpen (((cycloCharHom K n).ker : Subgroup K.absGal) : Set K.absGal) :=
  Subgroup.isOpen_mono (muFixer_le_ker_cycloCharHom K n) (isOpen_muFixer K (p ^ n))
```

★決め手は **`Subgroup.isOpen_mono`**（`Topology/Algebra/OpenSubgroup.lean:255`）——
「**開部分群を含む部分群は開**」。核が開であることを直接示す必要はなく、
既に開だと分かっている部分群（この木では `muFixer`）を核が含むことだけ言えばよい。

★★**仮説として置くときは、この形のほうが `Continuous` より弱い**（位相インスタンスを
要求しない）ので、後から本物の連続コホモロジーで具体化するときに困らない。

## #194 `MonoidHom.comap_ker` は逆向き —— `rw [← MonoidHom.comap_ker]`（2026-09-07、同上）

`(ψ.comp f).ker` を `ψ.ker.comap f` に書き換えたくて `rw [MonoidHom.comap_ker]` と書くと

```
Tactic `rewrite` failed: Did not find an occurrence of the pattern
  Subgroup.comap ?f (MonoidHom.ker ?g)
in the target expression
  IsOpen ↑(ψ.comp (↑α).toMonoidHom).ker
```

で落ちる。★mathlib の向きは `ψ.ker.comap f = (ψ.comp f).ker` なので **`←` が要る**。

```lean
theorem isOpen_ker_comp (K K' : PAdicLocalField p)
    (α : ContinuousMulEquiv K.absGal K'.absGal) (n : ℕ)
    (ψ : K'.absGal →* (ZMod (p ^ n))ˣ)
    (h : IsOpen ((ψ.ker : Subgroup K'.absGal) : Set K'.absGal)) :
    IsOpen (((ψ.comp (α : K.absGal ≃* K'.absGal).toMonoidHom).ker : Subgroup K.absGal)
      : Set K.absGal) := by
  rw [← MonoidHom.comap_ker, Subgroup.coe_comap]
  exact h.preimage α.continuous
```

★エラーメッセージが「探したパターン」を出してくれるので、
**そこに `comap` が出ていたら向きが逆**だと即断できる。

## #195 `Nat.card (groupCohomology.H2 A)` は**型が付く**が `Finite` は付かない —— 「書ける」を「使える」と読み違えない（2026-09-07、同上）

mathlib の離散群コホモロジーで

```lean
noncomputable example (A : Rep (ZMod m) G) : ℕ := Nat.card (groupCohomology.H2 A)   -- ★通る
example (A : Rep (ZMod m) G) [Finite G] [Finite A.V] :
    Finite (groupCohomology.H2 A) := by infer_instance                              -- ★落ちる
-- failed to synthesize instance of type class
--   Finite ↑(groupCohomology.H2 A)
```

★`groupCohomology.H2 A : ModuleCat k` には `CoeSort` があるので `Nat.card` は**書けてしまう**。
しかし `Finite` インスタンスは**有限群・有限係数でも無い**（2026-09-07 実測）。
すなわち `Nat.card (H2 A) = p ^ n` と書いても、`Nat.card` は無限のとき `0` を返すので
**主張が意図とずれたまま通ってしまう**危険がある。

★**教訓**: コホモロジーの「位数」を仮説や結論に書く前に、
**`Finite` インスタンスが引けるかを別行で `infer_instance` して確かめる**こと。
★同じ検査で分かったこと（在庫測定として再利用可）:
`groupCohomology.H1InfRes` は**次数 1 だけ**で `H2InfRes` は存在せず、
`groupCohomology.colimitIso`（有限群の塔から副有限群へ）も存在しない。
`e : G ≃* H` に沿った `H2 (Rep.res e.toMonoidHom A) ≅ H2 A` も `exact?` が閉じられない
（`groupCohomology.map` と `congr` から組めるはずだが束ねられていない）。

## #196 `rw [map_pow]` が `cyclotome` / `powTorsion` の境界で落ちる —— 片方の座標で `have` を書く（2026-09-07、Λ13 = Artin 写像の単数部分の同変性）

`cyclotome A m`（`Subgroup (TopologicalAbelianization A)`）と
`powTorsion (TopologicalAbelianization A) m` は**定義が同じ**（どちらも `(powMonoidHom m).ker`）だが、
`rw` は書き換え後に `instances` 透明度で型検査するので、両者をまたぐと落ちる。

```lean
refine (powTorsionCongr Θ (p ^ n)).injective ?_
rw [powTorsionCongr_abelianGalConj, map_pow]   -- ★落ちる
-- Tactic `rewrite` failed: Did not find an occurrence of the pattern ?f (?a ^ ?n)
-- Note: The target expression is not type-correct under the `instances` transparency level
```

★**直し方は 2 つある。**

(a) `map_pow` の**その場のインスタンス**を型を書いた `have` で作ってから `rw` する。

(b) ★**こちらが速い**: 等式そのものを**片方の座標だけ**で述べた補題に切り出し、
`(… ).trans (map_pow _ _ _).symm` のように**項で**繋ぐ（`exact` は既定透明度なので通る）。

```lean
exact (powTorsionCongr_abelianGalConj_eq_pow F n S hopen g E hE _).trans (map_pow _ _ _).symm
```

★同じ現象は `Subtype.val ⟨u, h⟩` と `u`（iota）でも起きる。
`rw` の連鎖の**最後に `rfl` を 1 行足す**だけで閉じることが多い。

## #197 `PowerSeries.map_map` と `Prod.fst_pow` は**無い**（2026-09-07、同上）

★測って分かったこと（`grep -n "PowerSeries.map_comp\|map_map" .cache/mathlib-index.txt`）:

| 書きたいもの | mathlib の実際 |
|---|---|
| `PowerSeries.map_map` | ★**無い**（`MvPowerSeries.map_map` だけ在る） |
| `Prod.fst_pow` | ★**無い**（`(a ^ k).1 = a.1 ^ k` は `rfl`） |

```lean
-- 合成は map_comp（RingHom の等式）を経由して rfl で降ろす
calc PowerSeries.map g (PowerSeries.map f x)
    = PowerSeries.map (g.comp f) x := by rw [PowerSeries.map_comp]; rfl
-- 直積の冪は rfl
have hfst : ∀ (a : Mˣ × N) (k : ℕ), (a ^ k).1 = a.1 ^ k := fun _ _ => rfl
```

## #198 `∀ a, α a = β a` 型の合同補題は `α` を `(α := …)` で明示する（2026-09-07、同上）

`topAbelianizationCME_congr {α β} (h : ∀ a, α a = β a)` のように
**関数適用の形でしか `α` が現れない**補題は、`exact` に渡すと `?m a =?= f (g a)` の
高階単一化になって落ちる。

```lean
exact topAbelianizationCME_congr (β := ContinuousMulEquiv.refl _)
  (fun a => (absGalFixedFieldCME F S hopen).symm_apply_apply a) w   -- ★落ちる
-- Type mismatch … but is expected to have type ?m.149 a = (ContinuousMulEquiv.refl …) a
```

★直し方: `(α := …)` も書く。そして `congr` 系は「`refl` に等しい」までしか言わないので、
**`topAbelianizationCME_refl` を `.trans` で継ぐ**のを忘れない。

```lean
exact (topAbelianizationCME_congr
  (α := (absGalFixedFieldCME F S hopen).trans (absGalFixedFieldCME F S hopen).symm)
  (β := ContinuousMulEquiv.refl _)
  (fun a => (absGalFixedFieldCME F S hopen).symm_apply_apply a) w).trans
  (topAbelianizationCME_refl w)
```

## #199 `IsWeierstrassFactorization.unique` は `refine … ?_` では単元部分の `?m` が決まらない（2026-09-07、同上）

`PowerSeries.IsWeierstrassFactorization.unique (H : g.IsWeierstrassFactorization f h) hg :
 f = g.weierstrassDistinguished hg ∧ h = g.weierstrassUnit hg` の `h`（単元部分）は
結論の `.1` からは決まらない。

```lean
refine (PowerSeries.IsWeierstrassFactorization.unique ?_ hg').1.symm   -- ★落ちる
-- don't know how to synthesize implicit argument `h`
```

★直し方: 分解を **3 引数すべて書いた `have`** にしてから `.unique` を当てる。

```lean
have H : (PowerSeries.map φ g).IsWeierstrassFactorization
    (Polynomial.map φ (g.weierstrassDistinguished hg))
    (PowerSeries.map φ (g.weierstrassUnit hg)) := ⟨…, …, …⟩
exact (H.unique hg').1.symm
```

★これで「Weierstrass 標準分解は係数のひねりと可換」（`map_weierstrassDistinguished`）が
**一意性 1 本**で出る。★`weierstrassDistinguished` の第 2 引数（`≠ 0` の証明）は
`subst` してから `rfl` で潰せる（`weierstrassDistinguished_congr`）。

## #200 `Rep k G` はもう `Action (ModuleCat k) G` ではない —— `Action.mkIso` が当たらない（2026-09-07、GroupCohomologyFinite）

`#print Rep` すると **フィールド `V` / `hV1` / `hV2` / `ρ` を持つ構造体**である。
したがって `Action.mkIso` は型が合わない（`M.V ≅ N.V` の `M N : Action V G` を要求する）。

★正しい作り方は `Rep.mkIso (e : ρ.Equiv σ) : Rep.of ρ ≅ Rep.of σ`。
`Rep.of A.ρ` は構造体 eta で `A` に defeq なので、
`Rep.res g (Rep.res f A) ≅ A` のような主張にそのまま使える。

```lean
noncomputable def resResIso {f : G →* H} {g : H →* G} (h : f.comp g = MonoidHom.id H)
    (A : Rep k H) : Rep.res g (Rep.res f A) ≅ A :=
  Rep.mkIso (Representation.Equiv.mk (LinearEquiv.refl k A) (fun x => by
    have hx : f (g x) = x := by rw [← MonoidHom.comp_apply, h]; rfl
    simp [hx]))
```

★`Representation.Equiv.mk (e : V ≃ₗ[A] W) (∀ g, ↑e ∘ₗ ρ g = σ g ∘ₗ ↑e)` の
絡み合い条件は `intro g` ではなく **`fun g => …`**（`∀ (g : G)` は既に剥がれている場合がある）。
実測で `intro` が `There are no additional binders` で落ちた。

## #201 `groupCohomology.congr` の `h ▸ φ` は使えない —— `f₁ f₂` を変数で量化して `subst` する（2026-09-07、同上）

mathlib の
`groupCohomology.congr (h : f₁ = f₂) (F) : F f₁ φ = F f₂ (h ▸ φ)`
は右辺に輸送 `h ▸ φ` が残るので、`φ` を `𝟙` と比較したいときに詰まる。

★直し方: **`f₁ f₂` を変数のまま量化した自前の合同則**を 1 本立てると `subst` が効く。
`φ` の比較は **台写像の pointwise 一致**だけでよい（`.hom` そのものは型が違うので比較できない
—— 実測で `Type mismatch: Rep.Hom.hom φ₂ has type (Rep.res f₂ A).ρ.IntertwiningMap B.ρ`）。

```lean
theorem map_congr {A : Rep k H} {B : Rep k G} {f₁ f₂ : G →* H} (h : f₁ = f₂)
    {φ₁ : Rep.res f₁ A ⟶ B} {φ₂ : Rep.res f₂ A ⟶ B} (hφ : ∀ x : A, φ₁.hom x = φ₂.hom x)
    (n : ℕ) : groupCohomology.map f₁ φ₁ n = groupCohomology.map f₂ φ₂ n := by
  subst h
  have : φ₁ = φ₂ := by ext x; exact hφ x
  rw [this]
```

★これで `e : G ≃* H` に沿った `groupCohomology A n ≅ groupCohomology (Rep.res e A) n` が
**全次数で** 0.5 秒で通る。

## #202 ★★`groupCohomology.map` 系の補題は暗黙引数のままだと `isDefEq` が heartbeat を焼く（2026-09-07、同上）

`rw [← groupCohomology.map_comp]` のあと `map_congr` を
`refine (map_congr (f₂ := …) (φ₂ := …) … ).trans (map_id n)` の形で当てると
**13 秒かけて `(deterministic) timeout at isDefEq`**。
`rw` で当てると今度は `φ₁` が勝手に単位射に unify されて
`Did not find an occurrence of the pattern`。

★直し方: **`A` `B` `f₁` `f₂` `φ₁` `φ₂` を全部明示**して `Eq.trans` で繋ぐ。
**0.5 秒**になる。#150（抽象核の穴は明示引数に）の `Rep` 版である。

```lean
exact Eq.trans (map_congr (A := Rep.res e.toMonoidHom A) (B := Rep.res e.toMonoidHom A)
    (f₁ := e.symm.toMonoidHom.comp e.toMonoidHom) (f₂ := MonoidHom.id G)
    (φ₁ := (Rep.resFunctor e.toMonoidHom).map (resMulEquivIso e A).hom ≫ 𝟙 _)
    (φ₂ := 𝟙 _) (by ext x; simp) (fun _ => rfl) n)
  (groupCohomology.map_id n)
```

★同じ理由で `intro x; rfl` は第 1 合成では通るが第 2 合成では
`timeout at whnf` になる。**明示引数にすれば両方 `rfl` で通る。**

## #203 `Finite (groupCohomology A n)` は 3 行で立つ —— mathlib に無いだけ（2026-09-07、同上）

`[Finite G]` `[Finite ↑A]` があっても `failed to synthesize Finite ↑(H2 A)`。
★中身は `Z^n ↪ C^n`（`iCocycles` は mono）と `Z^n ↠ H^n`（`π` は epi）だけである。

```lean
instance (A : Rep k G) (n : ℕ) [Finite G] [Finite A] : Finite ((inhomogeneousCochains A).X n) := by
  show Finite ((Fin n → G) → A); infer_instance
instance (A : Rep k G) (n : ℕ) [Finite G] [Finite A] : Finite (groupCohomology.cocycles A n) :=
  Finite.of_injective (ConcreteCategory.hom (groupCohomology.iCocycles A n))
    ((ModuleCat.mono_iff_injective _).1 inferInstance)
instance (A : Rep k G) (n : ℕ) [Finite G] [Finite A] : Finite (groupCohomology A n) :=
  Finite.of_surjective (ConcreteCategory.hom (groupCohomology.π A n))
    ((ModuleCat.epi_iff_surjective _).1
      (HomologicalComplex.instEpiHomologyπ (inhomogeneousCochains A) n))
```

★`Epi (groupCohomology.π A n)` は `infer_instance` では**出ない**が
`exact?` が `HomologicalComplex.instEpiHomologyπ` を即答する。
★`Finite ((inhomogeneousCochains A).X n)` も `infer_instance` では出ない。
**`show Finite ((Fin n → G) → A)` を 1 行挟むだけで出る**（`X n` が `ModuleCat.of` に
包まれているのでインスタンス探索が中に入れない）。

## #204 `ModuleCat` の同型から台の `Equiv` を作る —— `Iso.toEquiv` は `Type` 専用（2026-09-07、同上）

`CategoryTheory.Iso.toEquiv : {X Y : Type u} → (X ≅ Y) → X ≃ Y` なので
`X Y : ModuleCat k` には当たらない。`exact?` も 2.9 秒かけて閉じられない。

★直し方: **忘却関手で `Type` に落としてから** `toEquiv`。

```lean
theorem natCard_eq_of_moduleCatIso {k : Type u} [CommRing k] {X Y : ModuleCat.{v} k}
    (i : X ≅ Y) : Nat.card X = Nat.card Y :=
  Nat.card_congr (((forget (ModuleCat.{v} k)).mapIso i).toEquiv)
```

## #205 `rw` は `⇑↑Φ` と `⇑Φ` を別物として扱う —— 抽象核は `RingEquiv` 版を**別宣言**で用意する（2026-09-07、Λ14 = σ_g が Lubin-Tate 塔を運ぶ段）

抽象核を `(Φ : Ω →+* Ω)` で書き、消費側が `(Φ : Ω ≃+* Ω)` を持っていると、
消費側のゴールは `⇑Φ '' …`、補題の左辺は `⇑↑Φ '' …` になる。
**defeq なので `exact` は通るが `rw` は通らない**（`rw` は構文的一致しか見ない）:

```
Tactic `rewrite` failed: Did not find an occurrence of the pattern
  ⇑↑Φ '' ↑(IntermediateField.adjoin K.carrier ?S)
in the target expression
  ⇑Φ '' ↑(IntermediateField.adjoin K.carrier ↑(iteratedLubinTateTorsionPoints …))
```

★直し方: **`RingEquiv` 版のラッパを 1 本足す**。中身は `→+*` 版をそのまま渡すだけで、
項モード（`isDefEq`）なので通る。

```lean
theorem image_adjoin_of_semilinear_equiv {k Ω : Type*} [Field k] [Field Ω] [Algebra k Ω]
    (Φ : Ω ≃+* Ω) (σ : k ≃+* k)
    (hΦ : ∀ a : k, Φ (algebraMap k Ω a) = algebraMap k Ω (σ a)) (S : Set Ω) :
    Φ '' (IntermediateField.adjoin k S : Set Ω)
      = (IntermediateField.adjoin k (Φ '' S) : Set Ω) :=
  image_adjoin_of_semilinear (Φ : Ω →+* Ω)
    (image_range_algebraMap_of_semilinear (Φ : Ω →+* Ω) σ hΦ) S
```

★`show` で書き換えるより安い（`show` は目標全体を書き直すので長くなる）。

## #206 暗黙引数が未解決な項に `.injective` を付けると `function expected` になる（2026-09-07、同上）

`φ.injective (by rw [h, map_zero])` は `φ` の型が確定していれば通るが、
`(integerRingEquiv σ).injective (by …)` は通らない:

```
function expected
  RingEquiv.injective (integerRingEquiv σ)
```

原因: `integerRingEquiv : K.carrier ≃ₐ[ℚ_[p]] K'.carrier → 𝒪[K.carrier] ≃+* 𝒪[K'.carrier]`
の `K'` がまだメタ変数なので、`Function.Injective ⇑(integerRingEquiv σ)` が
`∀ ⦃a b⦄, … → …` へ**展開されない**（strict implicit の挿入が起きない）。

★直し方: **その場で書かず、既に在庫にある名前付き補題を呼ぶ**。
本ファイルでは `integerRingEquiv_ne_zero K σ hπne0`（`ArtinEquivarianceProof.lean`）が
そのまま当たった。★同じ理由で `maximalIdeal_eq_span_map (integerRingEquiv σ) hπmax` も
`maximalIdeal_eq_span_integerRingEquiv K σ hπmax` に置き換えると読みやすい
（どちらも **statement が同じなので proof irrelevance で下流と一致する**）。

## #207 ★★`@[simps]` の `_f` を `rw` すると「instances 透明度で型が合わない」（2026-09-07、InflationRestrictionH2）

`ShortComplex` を `def` で作り `@[simps X₁ X₂ X₃ f g]` を付けたあと、
`rw [H2InfRes_f] at hz` で `.f` だけを展開すると **hz が型不正になる**:

```
Note: The target expression is not type-correct under the `instances` transparency level
Application type mismatch: … has type
  … ⟶ groupCohomology (Rep.of A.ρ) 2
but is expected to have type
  … ⟶ (H2InfRes A S).X₂
```

`.X₂` は `def` のフィールドなので `instances` 透明度では開かない。
mathlib は `simp only [H1InfRes_X₁, H1InfRes_X₂, H1InfRes_f, …]` と**まとめて**書いて逃げているが、
そのあと `induction z using H2_induction_on` を打つと今度は

```
Dependent elimination failed: Failed to solve equation
  (H2InfRes A S).X₂.isAddCommGroup.toAddZeroClass.toZero.1 = …
```

で落ちる（`z` の型が `↥(H2InfRes A S).X₁` のままだから）。

★★直し方: **主定理を `ShortComplex` のフィールドではなく生の `map …` について述べる**。
`ShortComplex` 側へは `exact` で移す（`exact` は default 透明度なので通る）:

```lean
theorem inflation₂_injective … :
    Function.Injective (ConcreteCategory.hom
      (map (A := A.quotientToInvariants S) (B := A) (QuotientGroup.mk' S)
        (Rep.ofHom (A.ρ.quotientToInvariants_lift S)) 2)) := …

theorem mono_H2InfRes_f … : Mono (H2InfRes A S).f :=
  (ModuleCat.mono_iff_injective _).2 (inflation₂_injective A S h1)   -- ★これは通る
```

★同じ形で `(H2InfRes A S).Exact` も
`ShortComplex.moduleCat_exact_iff_ker_sub_range` → `intro c hc; exact <生の定理> c hc` で通る。

## #208 ★`Rep.ofHom` は暗黙引数を先に決めないと落ちる —— #202 の「timeout でない版」（2026-09-07、同上）

```lean
map (QuotientGroup.mk' S) (Rep.ofHom <| A.ρ.quotientToInvariants_lift S) 2
```
を**定理の statement に直接書く**と

```
Application type mismatch: The argument A.ρ.quotientToInvariants_lift S has type
  Representation.IntertwiningMap (MonoidHom.comp (A.ρ.quotientToInvariants S) (QuotientGroup.mk' S)) A.ρ
but is expected to have type
  Representation.IntertwiningMap (MonoidHom.comp (Rep.ρ ?m.36) (QuotientGroup.mk' S)) ?m.52
```

`structure` の中（`f := …`）では期待型が先に決まるので通るのに、statement では通らない。

★直し方: **名前付き暗黙引数を全部書く**。#202（`isDefEq` timeout）と同じ処方で、
症状だけが「timeout」でなく「型不一致」に変わったもの:

```lean
map (A := A.quotientToInvariants S) (B := A) (QuotientGroup.mk' S)
  (Rep.ofHom (A.ρ.quotientToInvariants_lift S)) 2
```

★おまけ: `Rep.ofHom` の codomain は `Rep.of A.ρ` であって `A` ではない。
`(B := A)` と書けば defeq で通る（default 透明度）が、`rw`（instances 透明度）では通らない。#207 参照。

## #209 抽象群の作用は `ρ : G → A →+ A` で受けると核が軽い（2026-09-07、同上）

コホモロジーの補題を「語彙ゼロの核」に落とすとき、作用の受け方は
`Representation k G V`（`k` が要る）でも `DistribMulAction`（instance を作る必要がある）でもなく
**`ρ : G → A →+ A` + 公理 2 本**が一番安い:

```lean
variable {G : Type u} {A : Type v} [Group G] [AddCommGroup A]
(ρ : G → A →+ A) (hone : ∀ a, ρ 1 a = a) (hmul : ∀ g h a, ρ (g * h) a = ρ g (ρ h a))
```

* `map_add` / `map_sub` / `map_zero` / `map_neg` が**ただで付く**（束ねた `A →+ A` だから）。
* 具体層からの代入は `fun g => (A.ρ g).toAddMonoidHom` の 1 行、
  しかも `(repAddHom A g) a = A.ρ g a` は **`rfl`**。
* `k` が消えるので universe も 2 つで済む。

★実測: この形にしたら `H²` の inflation-restriction の核 14 本が
**すべて 1 発（合計 0.3 秒）**で通った。`H¹(S,A) = 0` すら
`hH1 : ∀ f : G → A, (1-コサイクル条件) → ∃ a, ∀ s ∈ S, f s = ρ s a - a`
という**関数についての述語**に落ちる（`groupCohomology` が出てこない）。

## #210 「半線型な環同型がノルムを保つ」は最小多項式ではなく `spectralNorm_unique_field_norm_ext` で出す（2026-09-07、SemilinearRestriction）

`norm_algEquiv_eq`（木、`≃ₐ[K.carrier]` 専用）は「`σ x` は `x` と同じ最小多項式の根」で
証明されているので、**半線型（基礎体を集合としてしか保たない）には効かない**。
半線型版は mathlib の**一意性**に乗せると 20 行で出る:

```lean
-- f y := spectralNorm k Ω (Φ y) を AbsoluteValue に束ね、
spectralNorm_unique_field_norm_ext (f := pullbackSpectralAbs Φ)
  (fun c => by rw [hΦ, spectralNorm_extends, hσ]) y
```

材料は `spectralMulAlgNorm_def`（乗法性・劣加法性）と
`eq_zero_of_map_spectralNorm_eq_zero`（非退化）だけ。
★`Φ` の `k`-線型性をどこにも使わない。前提は
「`Φ` は `σ` について半線型」＋「`σ` は `k` のノルムを保つ」の 2 つ。
★★**系として「制限は等長 ⇒ 連続」が出るので、
`LinearMap.continuous_of_finiteDimensional`（線型性が要る）を使わずに済む。**
これで「cross-point instance bridging」（別々の `adjoinIntegers K x` /
`adjoinIntegers K (Φ x)` を橋渡しする段）が正面から越えられる。

## #211 `(fixedFieldAut S g).restrictScalars ℚ_[p]` を裸で書くと instance が落ちる（2026-09-07、同上）

`fixedFieldAut S g : ↥(IntermediateField.fixedField S) ≃ₐ[k] _` を
`.restrictScalars ℚ_[p]` した項を**引数位置に裸で置く**と、Lean は基礎型を
`↥(IntermediateField.fixedField S)` のまま取り、
`(fixedFieldLocalField F S hopen).carrier` 経由の instance
（`Algebra _ (fixedFieldLocalField F S hopen).closure` など）を見つけられず

```
failed to synthesize instance of type class
  Algebra (↥(IntermediateField.fixedField S)) (fixedFieldLocalField F S hopen).closure
```

で落ちる。★**直し方は名前付き引数で基礎型を先に固定する**:

```lean
conjSemilinearAlgEquiv (k := (fixedFieldLocalField F S hopen).carrier) ...
```

型注釈（`( ... : A ≃ₐ[ℚ_[p]] A)`）でもよいが、名前付き引数の方が短い。
★2 つの表示が defeq でも**インスタンス探索は syntactic な型で走る**、という #205 の親戚。

## #212 ★★`linear_combination (norm := abel)` の係数の符号は `2 •` が教えてくれる（2026-09-07、Transgression）

`linear_combination (norm := abel) e` が

```
⊢ 2 • X + (-2 • Y + (2 • Z + ...)) = 0
```

という**係数がすべて 2 の倍数**の残差で落ちたら、それは
「その項の符号だけが逆」という意味である。★残差は
`goal_diff - Σ cᵢ eᵢ` なので、正解が `+e` のときに `-e` を書くと
残差はちょうど `2 · (goal_diff)` になる。★**係数を反転すれば通る。**

実測（Transgression で 5 回）: `-h` → `h` / `e1 - e2` → `e2 - e1` /
`hys - E1 - E2` → `E1 - hys - E2`。★**残差を読んで手で解くより、
まず符号を反転して 1 往復（0.3 秒）で試す方が速い。**

★`linear_combination h` の規約は `hᵢ` について `(lhs - rhs)`。
`h : a = b + c` なら寄与は `a - b - c` である。★等式を `have` で作るとき
どちら向きに書いたかを忘れやすいのが原因。

## #213 ★`d₀₁` / `d₁₂` を当てた項は `abel` にとって別のアトム（2026-09-07、同上）

`hb : (d₀₁ A).hom b = ⇑x - ⇑y` から `congrFun hb s` を取ると

```
(ModuleCat.Hom.hom (d₀₁ (Rep.res S.subtype A))) b ⟨s, hs⟩ = ... ⟨s, hs⟩ - ... ⟨s, hs⟩
```

という形になり、`ρ s b - b` に**簡約されない**（defeq ではあるが構文が違う）。
`simp only [Pi.sub_apply]` を足しても `d₀₁ ... b ⟨s,hs⟩` は残るので、
`linear_combination`/`abel` は**それを 1 つのアトム**として扱い、閉じない。

★直し方は **`have` の型に欲しい形を書く**（defeq なので `congrFun` がそのまま通る）:

```lean
have hth : A.ρ s b - b = (x : G → A) s - f s := congrFun hb (⟨s, hs⟩ : S)
```

★これは `mapCocycles₁ f φ x ⟨s,hs⟩ = x ↑s` や
`⟨fun s => f ↑s, _⟩ ⟨s,hs⟩ = f s` にも効く（どれも `rfl`）。
★**「simp で潰す」より「型を書いて defeq に運ばせる」方が短い。**

## #214 `def` の中の `if g ∈ S then _ else _` は `Decidable (g ∈ S)` を要求する（2026-09-07、同上）

`Subgroup` の元判定は可判定でないので、`def` の中に書くと

```
failed to synthesize instance of type class Decidable (g ∈ S)
```

で落ちる（★`by classical` は**タクティク証明の中**でしか効かない）。

★`open scoped Classical` を足すより、**`if` を使わない式を探す方が良い**。
実例（正規化された左剰余類の代表元 `r 1 = 1`）:

```lean
noncomputable def normCosetRep (S : Subgroup G) (g : G) : G := cosetRep S g * (cosetRep S 1)⁻¹
```

★`if g ∈ S then 1 else cosetRep S g` と同じ 3 性質
（`r 1 = 1` / `(r g)⁻¹ g ∈ S` / `r (g s) = r g`）を満たすうえに、
★**可判定性も `S` の正規性も要らない**（`if` 版は正規性が要る場合がある）。

## #215 ★`LinearMap.range (ConcreteCategory.hom f)` は elaborate できない（2026-09-07、同上）

`f : X ⟶ Y`（`ModuleCat k`）について `LinearMap.range (ConcreteCategory.hom f)` と書くと

```
failed to synthesize instance of type class
  ConcreteCategory (Module k ↑X) (@LinearMap ...)
```

`Application type mismatch: ... has type Quiver.Hom.{0,1} X Y
 but is expected to have type Quiver.Hom.{0,0} X.isModule X.isModule`

になる。★`LinearMap.range` が期待する型から逆に `ConcreteCategory` を探しに行くため。

★直し方は **`ModuleCat.Hom.hom f` に書き換える**（`ConcreteCategory.hom f` と `rfl` で一致する）:

```lean
LinearMap.range (ModuleCat.Hom.hom (groupCohomology.map S.subtype (𝟙 _) 1))
```

★逆に「元に当てる」用途（`(ConcreteCategory.hom f) x`）はそのままでよい。
★**当てるときは `ConcreteCategory.hom`、`ker`/`range` を取るときは `ModuleCat.Hom.hom`。**

## #216 `rw [← map_add]` のあとの `congr 1` は goal を**閉じてしまう**（2026-09-07、同上）

`H1π ⟨f + f', _⟩ = H1π (⟨f,_⟩ + ⟨f',_⟩)` のような形は
`congr 1` が `Subtype` の値まで見て `rfl` で閉じる。
★そこに `exact Subtype.ext rfl` を続けると **`No goals to be solved`** になる。

★`congr 1` の後に何かを書く前に、**まず `congr 1` だけで通るか測る**（0.2 秒）。
★逆に閉じない場合だけ `apply Subtype.ext; funext …` を足す。

## #217 `induction … with | @H h` は同名の仮説を隠す（2026-09-07、同上）

```lean
theorem foo … (h : ∀ g j, P g j) … := by
  induction q1 using QuotientGroup.induction_on with | @H g =>
  induction q2 using QuotientGroup.induction_on with | @H h =>   -- ★ここで h が上書き
  rw [h g h]                                                     -- Function expected at h
```

★`QuotientGroup.induction_on` / `H1_induction_on` / `H2_induction_on` は
`with | @H x` で**名前を自分で決める**ので、仮説名とぶつかりやすい。
★仮説側を `hd` / `hFv` のように**前置き付き**にしておくと事故が減る。

## #218 ★★証明項を引数に取る作用を Subtype に包んだら、`rw` ではなく `exact`(定義的証明無関係)で閉じる（2026-09-07、ReciprocityLimitEquivariance）

```lean
-- psiRootAct U y := ⟨↑↑(unitActionQuotientLift … ↑y (mem_… y.2) (mem_adjoin_simple_self …) U), _⟩
theorem psiRootAct_mul … : psiRootAct a (psiRootAct b y) = psiRootAct (a * b) y := by
  apply Subtype.ext
  rw [coe_psiRootAct, coe_psiRootAct, coe_psiRootAct]   -- ★motive is not type correct
```

★理由: 内側の `psiRootAct b y` は評価点の位置にも、その点が `Λ_n` に属するという
証明の位置にも現れる。`coe_psiRootAct` で点だけを書き換えると証明の型が壊れるので
`rw` は motive を作れない(エラーも「e が a の性質に依存している」と言う)。

★★直し方: `rw` を全部やめて `exact` 1 本にする。

```lean
  apply Subtype.ext
  exact (lubinTateActionAtTorsionPoint_mul_eq_of_action … (v : 𝒪) (u : 𝒪) …).symm
```

★`Subtype` の `val` は `rfl` で潰れ、membership の証明は Prop の定義的無関係で
一致するので、`exact` の defeq 判定はそのまま通る(実測 1.4 秒)。
★同じ理由で `QuotientGroup.mk u * QuotientGroup.mk v` を `QuotientGroup.mk (u*v)` に
直す `rw [← QuotientGroup.mk_mul]` も要らない(`Quotient.map₂` なので `rfl`)。

★★一般則: 「証明項を引数に取る `def` を `Subtype` で包んだもの」の等式は
`simp`/`rw` ではなく `exact`(または `Subtype.ext` + `exact`)で閉じる。

## #219 ★木の Lubin-Tate 補題を `rw` に渡すときは明示引数を全部書く（2026-09-07、同上）

```lean
  rw [reciprocityUnits_absGalConjCME]     -- ★goal が 4 つ増える(hπmax / hf0 / hf1 / hf)
```

★理由: この木の Lubin-Tate 補題は `hπmax` `hπne0` `f` `hf0` `hf1` `hf` を明示引数で
持つ。左辺のパターンに現れないもの(`hf0` など `f` の性質)は `rw` がメタ変数のまま
残し、新しい goal として出てくる。
★見分け方: 残った goal が `⊢ (PowerSeries.coeff 0) f = 0` のように元の仮説そのものなら、これ。

★★直し方: `rw [lemma F S hopen hq g hπmax hπne0 f hf0 hf1 hf τ]` と全部書く。
★`rw` の後に `rfl` が自動で試されるが、`coe_fixedFieldIntegerAut` のような `rfl` 補題は
reducible 透明度では潰れないので、最後に `exact coe_… _` を足すこと。


## #220 ★`Module.DirectLimit.of R ι _ f i x` の族を `_` にすると instance が stuck（2026-09-07、CohomologyColimit）

```lean
    H2ColimToH2 A S hS (Module.DirectLimit.of k ι _ (h2Sys A S hS) i x) = inflH2Lin A (S i) x
-- error: typeclass instance problem is stuck
--   Preorder (ModuleCat k)
```

★理由: `Module.DirectLimit.of R ι G f i` の `G : ι → Type*` を `_` にすると、
Lean は `f` の型から `ι` を決める前に `G` を `ModuleCat k` そのものと読もうとして、
添字の `Preorder` を `ModuleCat k` に探しに行く。

★★直し方: **族を全部書く**。

```lean
    Module.DirectLimit.of k ι
      (fun i => (groupCohomology (A.quotientToInvariants (S i)) 2 : ModuleCat k))
      (h2Sys A S hS) i x
```

★同じ形で `Module.DirectLimit.lift_of _ _ x` も落ちる。`lift_of g Hg x` と
`g`・`Hg` を明示すれば通る(#202 の系)。

## #221 `Module.DirectLimit` は `DirectedSystem` を要求しない（2026-09-07、同上）

★実測: `Module.DirectLimit G f` の定義にも `Module.DirectLimit.exists_of` にも
`[DirectedSystem …]` は要らない。要るのは

* 型を作るのに `[DecidableEq ι]`（★**型が instance に依存する**ので、
  `Classical.decEq` を `letI` で差し込まず**引数で受ける**こと）、
* `exists_of` / `induction_on` に `[Nonempty ι]` と `[IsDirectedOrder ι]`。

★したがって「有向系であること(`map_self` / `map_map`)」を証明しなくても colimit の
定理は書ける。★それでも別途証明しておくと、系が本当に関手的であることの自己検査になる
（本ファイルは `resH1Step_self` / `resH1Step_trans` / `inflH2Step_self` /
`inflH2Step_trans` を置いた。どれも `congr 1` + `QuotientGroup.induction_on` + `rfl`）。

## #222 `t : ↥(⊥ : Subgroup G)` に `simpa using t.2` を当てると仮説が `True` に潰れる（2026-09-07、同上）

```lean
  have h1 : t = 1 := Subtype.ext (by simpa using t.2)
-- error: term t.property has type True but is expected to have type t = 1
```

★理由: simp は `↑t ∈ ⊥` を `↑t = 1` にしたあと、**`t.2` 自身で閉じて `True` にする**。

★★直し方: simp を使わず `Subgroup.mem_bot.1 t.2` と書く。

```lean
  have h1 : t = 1 := Subtype.ext (Subgroup.mem_bot.1 t.2)
```

★ただし `ht : (t : G) ∈ S j` のように**別の仮説**を書き換えるときは `simpa using ht` で通る。
違いは「その仮説が結論そのものになっているか」である。

## #223 `coe_mapCocycles₁` のあとの `Rep.Hom.hom (𝟙 _)` は simp で消えない（2026-09-07、同上）

```lean
  ext t
  rw [coe_mapCocycles₁]
  simpa using hz t          -- ★NG: goal に (Rep.Hom.hom (𝟙 (Rep.res …))) が残る
```

★★直し方: `show` で恒等射を**書き下してから** `rw` する。

```lean
  ext t
  rw [coe_mapCocycles₁]
  show (Rep.Hom.hom (𝟙 (Rep.res (Subgroup.inclusion hle) (Rep.res S.subtype A))))
      ((z : ↥S → A) (Subgroup.inclusion hle t)) = (0 : ↥T → A) t
  rw [hz (Subgroup.inclusion hle t) t.2]
  simp
```

★一方 `ConcreteCategory.hom (mapCocycles₁ f (𝟙 _)) z` を**関数として**述べた補題
（`coe_mapCocycles₁_inclusion`）を作っておくと、`rw` 1 発で `rfl` まで閉じる。
★★**恒等射が絡む `mapCocycles` は「値の補題」を 1 本先に作る**のが安い。


## #224 ★★★明示引数が残っている補題に `.mp` / `.1` を打つ罠 —— ★**1 つの罠に 2 つの顔がある**（2026-09-07、メタ第 26 回の実測）

★★**この節は「逆引き」専用である。**★下の 2 行は **Lean の出力を逐語で引いたもの**。

```
Invalid projection: Projections cannot be used on functions, and <補題名> ?m.216 has function type
Invalid field `mp`: The environment does not contain `Function.mp`, so it is not possible to project
```

★★**どちらも同じ原因**: `Iff` を返す補題に **明示引数がまだ残っている**のに
`.mp` / `.mpr` / `.1` / `.2` を打った。★項はまだ `Iff` でなく**関数**なので、
Lean は①射影を拒む（前者）か②`Function` 名前空間の場を探しに行って失敗する（後者）。
★**`?m.NNN`（メタ変数）が入っていたらこの罠である**（実測で 40 件中 27 件）。

★**直し方**: ★★**明示引数を全部書く**。`_` で済ませない。
★探す道具は `.cache/mathlib-index.txt`（statement ごと）か `lean_check` の `#check`。

★★**実測（なぜこの節を作ったか）**: ★合わせて **40 件**、
★**08-27 から 09-07 まで毎日**、★**09-07 だけで 6 件**。
★同じ現場が両方の顔を出している（`mem_fixingSubgroup_iff` は 09-06 に後者、09-07 に前者）。
☆★**この罠は既に #72 / #104 / #114 / #125 の 4 節に散っていたが、
4 つとも逐語でなかった**（引用符が `'`、途中が `...`）。
⇒ ★**人の grep も機械の照合も当たらなかった。それで再発し続けた。**

★★★**作法（ここが一般化できる）**: ★**1 つの罠が複数のエラー文を持つときは、
その全部を逐語で並べること。**★言い換えや省略記号を入れた瞬間に索引から落ちる。

---

## #225 ★★`Rep k G` のまま次数 `n` の `d_comp_d` を使うと**宇宙が合わない**（2026-09-07、ContinuousCochain）

```
Application type mismatch: The argument
  A
has type
  Rep.{u_1, 0, 0} k G
but is expected to have type
  Rep.{u_1, u_1, u_1} ?m.33 ?m.34
in the application
  groupCohomology.inhomogeneousCochains.d_comp_d n A
```

★**原因**: `inhomogeneousCochains.d`（★**根の名前空間**。`groupCohomology.` が付かない）は
係数の宇宙が自由（`Rep.{v, u, u}` を受ける）だが、
`groupCohomology.inhomogeneousCochains` と `groupCohomology.inhomogeneousCochains.d_comp_d` は
`Rep.{u,u,u}` を要求する。`variable {k G : Type}` + `(A : Rep k G)` と書くと、
その宣言で `v` が **auto-bound の `u_1` のまま**になり、上のエラーになる。

★**直し方**: ★**`(A : Rep.{0} k G)` と書く**。★**ファイル内の全宣言で統一すること**
（片方だけ直すと今度は宣言同士で型が合わなくなる）。

★**なぜ低次では出ないか**: statement に `groupCohomology A n` や `cocycles₁ A` が
現れる宣言では `v = 0` が**強制される**ので起きない。
⇒ ★**「低次の補題では一度も出ず、次数 `n` の複体を触った瞬間に出る」**。

## #226 ★★`ModuleCat` 越しの `f = 0` を `congrFun` で開くと `Pi.zero_apply` が効かない（2026-09-07、同上）

```
Type mismatch: After simplification, term
  h
 has type
  (A.ρ g) (f default) + -f default = 0 fun x => g
but is expected to have type
  (A.ρ g) (f default) = f default
```

★`0 fun x => g` が消えない。`simp` は
`This simp argument is unused: Pi.zero_apply` と言って**当たらなかったことを申告する**。
原因は `0` が `↑(ModuleCat.of k ((Fin 1 → G) → ↑A))` の zero であって、
`Pi.zero_apply` の左辺（`OfNat.ofNat 0 x` の Pi 版）と**構文的に一致しない**こと。

★**直し方**: `congrFun` の結果に**Pi 型で型を付け直す**。

```lean
have h : ModuleCat.Hom.hom (inhomogeneousCochains.d A 0) f (fun _ => g)
    = (0 : (Fin (0 + 1) → G) → ↥A) (fun _ => g) := congrFun hker _
```

★goal 側（`d f = 0` を示す方）は先に
`show ModuleCat.Hom.hom (inhomogeneousCochains.d A 0) f = (0 : (Fin (0 + 1) → G) → ↥A)`
を打ってから `funext` する。

## #227 ★`fun g => Φ g` の形の目標は `rw` が当たらない（beta 未簡約）（2026-09-07、同上）

```
Tactic `rewrite` failed: Did not find an occurrence of the pattern
  j.contractNth ?m.206 g i
in the target expression
  ((fun g => j.contractNth (fun x1 x2 => x1 * x2) g) g i)⁻¹ *
      (fun g => j.contractNth (fun x1 x2 => x1 * x2) g) (fun i => g i * h i) i ∈
    V
```

★述語の定義に `Φ : (κ → G) → (ι → G)` を**ラムダのまま**代入すると、
目標に `(fun g => …) g` が残り、`rw` の照合が落ちる（`simp only []` も進まないことがある）。

★**直し方**: `rw` の前に **`show` で beta 簡約した形を書く**。

```lean
show (Fin.contractNth j (· * ·) g i)⁻¹ *
  Fin.contractNth j (· * ·) (fun i => g i * h i) i ∈ V
```

★**同じ日にもう 1 つの顔**: `QuotientGroup.eq_one_iff.2 (hh i)` が
`Invalid projection: Projections cannot be used on functions, and QuotientGroup.eq_one_iff has function type ∀ (x : ?m.153), ↑x = 1 ↔ x ∈ ?m.155`
で落ちた（★**#224 の再発。登録後に 1 件**）。直し方は `(QuotientGroup.eq_one_iff (h i)).2 (hh i)`。

## #228 ★`convert` が生む「位相インスタンスの等式」ゴールは**関数のゴールより先**に出る（2026-09-07、ReciprocityDatumIndependence）

```
⊢ instTopologicalSpaceSubtype = instUniformSpaceSubtype.toTopologicalSpace
```

★`HasSum` を `Subtype`（`ValuationSubring` の coe）の上で `convert h using 2 with d` すると、
**関数の等式**のほかに**位相インスタンスの等式**が goal に出る。`rfl` で閉じるが、
★**出る順序が「インスタンス → 関数」**なので、bullet を書かずに `rw` を続けると
`rw` がインスタンスの goal に当たって
`Did not find an occurrence of the pattern … in the target expression instTopologicalSpaceSubtype = …`
で落ちる。

★**直し方**: bullet を 2 つ書き、**先に `· rfl`**、後に本体。

```lean
convert hs3 using 2 with d
· rfl
· rw [PowerSeries.coeff_map, …]
```

## #229 ★`PowerSeries.HasEval` の**評価点を自分で書かない**（2026-09-07、同上）

`hasEval_coe_of_norm_lt_one K hlam : HasEval ⟨↑lam, ⋯⟩` の `⋯` は
`(mem_closureCompletionInt K _).mpr (by …)` という**書き下しにくい項**である。
`(hasEval_… ).elim (fun _ _ => trivial)` のような誤魔化しは
`Invalid field elim` で落ちる（`HasEval` は `Filter` の含意であって構造体ではない）。

★**直し方**: 点を書かず、`refine` で**補題の型から決めさせる**。

```lean
have h1 : PowerSeries.aeval hv θ = evalAt K θ hμnorm := by
  refine aeval_congr_point K hv (hasEval_coe_of_norm_lt_one K hμnorm) ?_ θ
  apply Subtype.ext
  exact …
```

## #230 ★係数環が違う 2 つの `aeval` を繋ぐには `hasSum_aeval` + `HasSum.map`（2026-09-07、同上）

`PowerSeries.comp_aeval (ha) (hε) : ε.comp (aeval ha) = aeval (ha.map hε)` は
**係数環 `R` が両辺で同じ**ときにしか使えない。
`aeval (R := 𝒪_K)`（`adjoinIntegers K x` への評価）と
`aeval (R := 𝒪_{K̂^ur})`（`𝒪_{ℂ_K}` への評価、`map (baseIntHom K) F` を入れる）を
繋ぎたいときは `comp_aeval` が当たらない。

★**直し方**: `PowerSeries.hasSum_aeval` → `HasSum.map φ hφ` → `HasSum.unique`。
係数の差は `(baseIntHom K a) • w = a • w`（`IsScalarTower.algebraMap_smul` の後 `rfl`）で吸収する。
★**これで `ContinuousSMul 𝒪_K 𝒪_{ℂ_K}` を組まずに済む**（`aeval (R := 𝒪_K)` を
`𝒪_{ℂ_K}` の上で一度も書かない）。

## #231 ★`eλ` は識別子にできない（2026-09-07、同上）

```
error: unexpected token 'λ'; expected ')'
```

`λ` は Lean 4 の予約トークンなので、`eλ`・`heλ` のような**λ を含む識別子は作れない**
（`θ`・`π`・`ϖ`・`σ` は使える）。★`elam`・`helam` に直す。

## #232 ★★★`grep -i` の「部分文字列に飲み込まれる」形 —— 正規底定理を「mathlib に無い」と誤判定した（2026-09-08、pGC Prop 2.1）

★**Lean のエラーではなく測定の失敗形**だが、実害は同じ（`.absent` を 1 件でっち上げるところだった）。

やったこと:

```
grep -in "normalBasis\|normal_basis" .cache/mathlib-index.txt | head -30
```

返ってきたのは `OrthonormalBasis`・`orthonormalBasis` **ばかり 30 行**で、
`Orthonormal**Basis**` の中に `normalBasis` が部分文字列として入っているため
本命が `head` の外に押し出されていた。★実際には在る:

```
Mathlib/FieldTheory/Galois/NormalBasis.lean
  IsGalois.normalBasis (K L) : Module.Basis Gal(L/K) K L
  IsGalois.normalBasis_apply (e : Gal(L/K)) : normalBasis K L e = e (normalBasis K L 1)
```

★**直し方（3 手、どれも 1 秒）**:

1. **大文字小文字を区別して型名で引く**: `grep -n "NormalBasis" .cache/mathlib-index.txt`
2. **飲み込む語を除く**: `... | grep -v -i orthonormal`
3. **`head` を付けるなら `wc -l` も見る**（30 行で切った時点で「全部見た」と思わない）

★同型の飲み込み例: `Basis`↔`OrthonormalBasis`/`HilbertBasis`、`Norm`↔`Enorm`/`Seminorm`、
`Inv`↔`Invariant`/`Involutive`、`Gal`↔`Galois`/`generalized`。

★★見つけたあとの最初のエラーは**不在ではなく import 漏れ**（#68）だった:

```
error(lean.unknownIdentifier): Unknown constant `IsGalois.normalBasis`
```

`import Mathlib.FieldTheory.Galois.NormalBasis` の 1 行で消えた。

## #233 `DistribMulAction` の instance を `where` の中で直接 `simp` すると `•` が開かない（2026-09-08、同上）

```
error: `simp` made no progress
error: unsolved goals
⊢ (g • 0) x = 0
error: unsolved goals
⊢ (g • (f₁ + f₂)) x = (g • f₁) x + (g • f₂) x
```

定義中の instance の `smul` フィールドは、その instance 自身の証明フィールドからは
**`simp` が展開できない**（インスタンスがまだ環境に無い）。

★**直し方**: 作用を先に `def` で切り出し、`smul := translate` と書いてから
各フィールドを `show` で開く。

```lean
def translate (g : G) (f : LocallyConstant G M) : LocallyConstant G M :=
  LocallyConstant.comap (leftTranslation g) f

scoped instance : DistribMulAction G (LocallyConstant G M) where
  smul := translate
  one_smul f := by ext x; show f (1⁻¹ * x) = f x; simp
  smul_add g f₁ f₂ := by
    ext x; show (f₁ + f₂) (g⁻¹ * x) = f₁ (g⁻¹ * x) + f₂ (g⁻¹ * x); simp
```

## #234 `rw [h]`（`h : α = α'`）は依存する項があると motive エラー（2026-09-08、同上）

```
Tactic `rewrite` failed: motive is not type correct:
  fun _a => _a g • e.toAddEquiv a = α' g • e.toAddEquiv a
Error: Application type mismatch: The argument
  e
has type
  SemilinearAddEquiv α A B
but is expected to have type
  SemilinearAddEquiv _a A B
```

★**直し方**: `rw` ではなく `subst`。

```lean
def ofEq {α α' : G ≃* G'} (h : α = α') (e : SemilinearAddEquiv α A B) :
    SemilinearAddEquiv α' A B := by
  subst h; exact e
```

## #235 `simpa using h` は `e.symm.toAddEquiv` と `e.toAddEquiv.symm` を同一視しない（2026-09-08、同上）

```
error: Type mismatch: After simplification, term
  h2
 has type
  @Eq B (eB.symm.toAddEquiv (α g • e.toAddEquiv (eA.toAddEquiv a)))
    (α g • eB.symm.toAddEquiv (e.toAddEquiv (eA.toAddEquiv a)))
but is expected to have type
  @Eq B (eB.toAddEquiv.symm (α g • e.toAddEquiv (eA.toAddEquiv a)))
    (α g • eB.toAddEquiv.symm (e.toAddEquiv (eA.toAddEquiv a)))
```

`symm` を自分で定義した構造では `e.symm.toAddEquiv` と `e.toAddEquiv.symm` は
**defeq だが構文的に別**。`simpa` は構文で照合するので落ちる。

★**直し方**: `have` に**明示的な型を書いて** defeq で受け、`exact` する。

```lean
have h2 : eB.toAddEquiv.symm (α g • c) = α g • eB.toAddEquiv.symm c :=
  eB.symm.map_smul (α g) c
```

## #236 ★★MCP の基準環境は**並行セッションと共有**される —— 自分の `lean_start` が上書きされる（2026-09-08、同上）

`lean_start` が「起動して import を読み込んだ (107.7 秒)」と**成功を返した直後**に、

```
error: failed to synthesize instance of type class
  AddCommGroup (LocallyConstant X Z)
error: Unknown constant `LocallyConstant.comapMonoidHom`
```

が出た。`Mathlib.Topology.LocallyConstant.Algebra` は確かに import 指定していたのに、である。
`lean_status` を見ると:

```
imports: ABC3.Found.PGC.AbsGalRamificationFiltration, ABC3.Found.PGC.LubinTateClosureTopology,
         ABC3.Found.PGC.LubinTateReciprocityMapLimitSurjective
建て直し 1 回
```

★**自分が指定した 3 本が 1 本も入っていない**（別セッションの import で建て直されていた）。

★**直し方**: 「Unknown constant なのに import したはず」と思ったら
**まず `lean_status` で imports を見る**。違っていたら**再起動を試みず**
`node tools/leanfile.mjs lean/ABC3/.../Foo.lean` へ切り替える（11〜13 秒/往復、olean を書かない）。

## #237 `Basis.equivFun` の往復を `simp` で潰すと総和に展開されて閉じない（2026-09-08、同上）

```
error: unsolved goals
⊢ ∑ x, f x * (Finsupp.single x 1) ρ = f ρ
```

`b.equivFun (b.equivFun.symm f) ρ = f ρ` に `simp` を当てると、
`Basis.equivFun_symm_apply`（`∑ i, x i • b i`）と `repr_self` が先に発火して
**総和に展開されてから**戻れなくなる。

★**直し方**: `rw [LinearEquiv.apply_symm_apply]`（往復をまず消す）。

## #238 `IsCompact.nonempty_iInter_of_directed_nonempty_isCompact_isClosed` は `ι : Type*` しか取らない（2026-09-08、pGC §2 分岐入力）

添字を `{ι : Sort*}` で書くと、**最後の引数だけ**が通らない:

```
error: Application type mismatch: The argument
  htcl
has type
  ∀ (i : ι), IsClosed (t i)
but is expected to have type
  ∀ (i : ?m.281), IsClosed (?m.284 i)
in the application
  IsCompact.nonempty_iInter_of_directed_nonempty_isCompact_isClosed ?m.284 ?m.285 ?m.287 (fun i => ?m.294) htcl
```

★`?m.284 : ?m.281 → Set Γ` の `?m.281` が `Type` 宇宙に固定されているので `Sort*` の `ι` が入らない。
**直し方: `{ι : Type*}` にする**（`Sort*` にする理由は無い）。
★「4 つ目までは通るのに 5 つ目で落ちる」ので、自分の `htcl` の側を疑って時間を溶かしやすい。

## #239 派生した族を渡すと、ユニフィケーションが「族そのもの」を間違える（2026-09-08、同上）

`{T : ι → Subgroup Γ}` を持つ補題を `F.S (T i) v` の形で使うと:

```
error: Application type mismatch: The argument
  hdir
has type
  Directed (fun x1 x2 => x1 ≥ x2) T
but is expected to have type
  Directed (fun x1 x2 => x1 ≥ x2) fun i => F.S (T i) v
```

★最後の引数（`↑S * ↑(T i) = …`）から `?T := fun i => F.S (T i) v` と**先に**拉致されている。
**直し方: `(T := T)` を名前付きで渡す**。★直したあとは
`IsClosed ↑(F.S (T i) v)` を渡していた所が `IsClosed ↑(T i)` に変わるので、

```
error: Type mismatch
  isClosed_S F (hbase i) v
has type
  IsClosed ↑(F.S (T i) v)
but is expected to have type
  IsClosed ↑(T i)
```

が続けて出る。**2 段階で直る**（`Subgroup.isClosed_of_isOpen _ (F.isOpen_base (hbase i))`）。

## #240 `Nat.ceil_eq_iff` の右辺は切り捨て引き算 `↑(n - 1)`（2026-09-08、同上）

```
error: Application type mismatch: The argument
  h1
has type
  ↑n - 1 < v
but is expected to have type
  ↑(n - 1) < v
in the application
  And.intro h1
```

★原典の丸め「`n − 1 < v ≤ n`」を `Nat.ceil_eq_iff` に渡すとここで落ちる。
**直し方: `1 ≤ n` を出して `rwa [Nat.cast_sub hn, Nat.cast_one]`**。
★`n = 0` は `Nat.ceil_eq_iff` の側条件（`n ≠ 0`）で弾かれるので、
**`rcases Nat.eq_zero_or_pos n` で先に場合分けする**。`v ≥ 0` なら `n = 0` のとき
`v = 0` なので `⌈v⌉₊ = 0` で合う。

## #241 ★★★Python の `io.open(p,'w')` は**書く前に空にする** —— `UnicodeEncodeError` でファイルが消える（2026-09-08、同上）

★**Lean のエラーではなく、木を壊した実害の形**である。`tools/lean-idioms.md`（10468 行）を
**0 バイトにした。**

```
Traceback (most recent call last):
  File "...add_idioms.py", line 97, in <module>
    io.open(p, 'w', encoding='utf-8', newline='\n').write(s)
UnicodeEncodeError: 'utf-8' codec can't encode characters in position 236660-236661: surrogates not allowed
```

原因は、Write ツールに渡した Python の中に `𝒪`(U+1D4AA) を
**`\ud835\udcaa` というサロゲート対のエスケープ**で書いたこと。Python 3 の文字列では
これは「対」にならず `'\ud835'` 単体（lone surrogate）になり、`utf-8` で書けない。
★**`'w'` はファイルを開いた瞬間に truncate する**ので、`.write` が例外を投げた時点で
中身は失われている。

★**直し方（3 つとも守る）**:
1. ★**Markdown/Lean への追記は Edit ツールでやる**（Edit は失敗すれば何も書かない）。
2. Python でやるなら **BMP 外の文字はエスケープせず直に書く**（`𝒪`・`ℝ`）。
   どうしてもエスケープするなら `\U0001D4AA`（大文字 `U` + 8 桁）。
3. Python でやるなら **一時ファイルに書いてから `os.replace`**（原子的に差し替える）。

★★**復旧の道（実際にこれで全部戻した）**: 未 commit の変更でも、
書き込んだ agent の**サブエージェント・ログ**に `Edit` の `old_string` / `new_string` が
そのまま残っている。

```
C:\Users\Aruta\.claude\projects\D--Math-ABC3\<session>\subagents\agent-*.jsonl
grep -l '## #232' *.jsonl
```

で当たりを引き、`tool_use` の `input` を JSON として取り出して同じ Edit を再適用すればよい。
★`git checkout` は **HEAD に戻すだけ**なので、未 commit の他 agent の追記を失う。先にログを見ること。

## #242 主張の中の `letI := …hfd` は、インスタンスが既に合成できると**消える**（2026-09-08、Yoshida Prop 6.14 の鋭い形）

構造体のフィールド（`(psiGenSeq … m).hfd : FiniteDimensional …`）を主張に持ち込もうとして

```lean
example … : letI := (psiGenSeq K … m).hfd
    reciprocityMap K … (psiGenSeq K … m).pt … σ = … := by
  intro _
  …
```

と書くと：

```
Tactic `introN` failed: There are no additional binders or `let` bindings in the goal to introduce
```

★原因は「`letI` が効かなかった」ではなく **`letI` が要らなかった**。
インスタンス探索が `FiniteDimensional K.carrier K⟮(psiGenSeq … m).pt⟯` を自力で見つけたので
`letI` が elaboration で潰れ、`intro` する束縛子が無くなっている。

★**直し方**: `letI` も `intro` も書かない。★**先に測ること**——

```lean
example … : True := by
  have h1 : FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier
      ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).pt} : Set K.closure)) := by infer_instance
  trivial
```

が 0.26 秒で通れば、主張に `[FiniteDimensional …]` 束縛子を足す必要も無い。
★この木では `Fintype (K⟮x⟯ ≃ₐ[K.carrier] K⟮x⟯)` も同様に合成できる。

## #243 `show` の中の `_` は「明示引数の証明」を埋められない（2026-09-08、同上）

`MulEquiv.ofBijective` の外側を剥がすために

```lean
    show galoisUnitReciprocityMap K … (M + 1) (by omega) _ _ _ _
        (algEquivRestrictSelf K … (M + 1) (by omega) _ _ _ _ σ) = _
```

と書くと：

```
don't know how to synthesize placeholder for argument `hxψ`
context:
⊢ (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf M).pt ∈
    iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (M + 1) ⋯
```

★`show` は目標との**構文的**な照合をしないので、`_` に入る `hxψ : x ∈ …`（Prop だが明示引数）を
逆算できない。★**全部書き下す**か、次の形にする（こちらが安い）：

```lean
theorem galoisUnitReciprocityEquiv_apply (y : …) :
    galoisUnitReciprocityEquiv K … n hn x hxψ hxn hmem y
      = galoisUnitReciprocityMap K … n hn x hxψ hxn hmem y := rfl
```

を `variable` 節（`(n) (hn) (x) (hxψ) (hxn) (hmem)`）の中で 1 本立てて `rw` する。
★引数が 13 個あるとき、`show` を 1 回書くより `_apply` 補題 1 本のほうが短い。

## #244 `apply e.injective` のあとの `rw` は「defeq だが構文が違う」ところで止まる（2026-09-08、同上）

`reciprocityMapLimitFamily` は `match n with | 0 => … | n+1 => principalUnitsQuotientEquiv … (reciprocityMap …)`
という `def` である。`apply (principalUnitsQuotientEquiv …).injective` して
`rw [principalUnitsQuotientEquiv_apply_mk, …, hkey]` まで進めると：

```
unsolved goals
⊢ (principalUnitsQuotientEquiv K hπmax (m + 1) ⋯)
      (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (m + 1) ⋯ (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt ⋯ ⋯ ⋯ σ) =
    reciprocityMapLimitFamily K hq hπmax hπne0 f hf0 hf1 hf σ (m + 1)
```

★両辺は **defeq**（右辺を `m+1` で展開すると左辺そのもの）だが、`rw` は最後に
`Eq.refl` の**構文的**照合しかしないので閉じない。★**`rfl` を 1 行足すだけ**でよい
（`rfl` は default transparency なので `def` の `match` を展開する）。

## #249 `.cache/mathlib-index.txt` は**取りこぼす** —— 索引に無いことは不在の証拠にならない(2026-09-08)

☆★**本体の追記(2026-09-08)**: ★**原因を特定して直した。**
mathlib が新しいモジュールシステムの `public` 修飾子を使い始めており
(`public noncomputable def normalBasis` / `public theorem normalBasis_apply`)、
★`tools/decl-index.mjs:58` の `MODS` に `public` が無かったため
★**`public` の付いた宣言が丸ごと索引から落ちていた。**
★`public` と `nonrec` を足して作り直した: ★**248,636 → 249,481 宣言(+845)**、
`IsGalois.normalBasis` が引けるようになった。
★**それでも「索引に無い ⇒ 不在」は言えない。**`.absent` を書く前に `#check @Foo` で確かめること。
☆★**この節は元々 `## 489.` と書かれており、`^## #` に当たらないので
`idiom-recur.mjs` から不可視だった。**★**番号は既存の最大＋1 にすること。**

**現象**: `IsGalois.normalBasis`(有限次 Galois 拡大の正規底、
`Mathlib/FieldTheory/Galois/NormalBasis.lean`)を索引で引くと**出ない**。

```
$ grep -n "normalBasis" .cache/mathlib-index.txt | head
16238:def  Complex.isometryOfOrthonormal  Analysis/InnerProductSpace/PiL2.lean:902 ...
...(OrthonormalBasis ばかり。IsGalois.normalBasis は 1 件も無い)
$ grep -n "Galois/NormalBasis.lean" .cache/mathlib-index.txt
236152:theorem exists_linearIndependent_algEquiv_apply_of_finite ...
236153:theorem exists_linearIndependent_algEquiv_apply_of_infinite ...
```

★**ファイル自体は索引されているのに、その中の `noncomputable def normalBasis` だけが
落ちている。**「同ファイルから 2 件出るのだから網羅されているはず」という推論は誤り。

**確かめ方(11 秒)**: スクラッチに 2 行書いて `node tools/leanfile.mjs` に投げる。

```lean
import Mathlib.FieldTheory.Galois.NormalBasis
#check @IsGalois.normalBasis
```
```
IsGalois.normalBasis : (K : Type u_1) → (L : Type u_2) → [inst : Field K] → … → Module.Basis Gal(L/K) K L
```

**How to apply**: 「mathlib に無い」と書く前に、★索引の grep だけで済ませない。
(i) 名前空間で grep、(ii) `#check @<推定名>` を 1 回投げる、(iii) `exact?`。
★2026-09-07 に Y21 が 4 件を「不在」と誤報告した件(#117(ii))の**逆側の失敗形**である
——あちらは grep の書き方、こちらは**索引そのものの欠落**。

★同じ 2026-09-08 の測定で、**本当に不在**だったものも記録しておく:
Krasner の補題は**在る**(`IsKrasner` / `IsKrasner.krasner`,
`Mathlib/Analysis/Normed/Field/Krasner.lean`)が、その系
「局所体の与えられた次数の拡大は有限個」「有限次部分拡大は可算」は
`grep -n "finite_extensions\|countable.*IntermediateField" .cache/mathlib-index.txt` が
**0 件**で、本木にも無い。★「在る/無い」は 1 語ではなく**主張の形**で測ること。

## #245 ★★余計な `haveI : Fintype … := Fintype.ofFinite _` を書くと、**見た目が同一の項に `rw` が当たらない**（2026-09-08、pGC「Art(Γ^n)=U^n」）

```
error: Tactic `rewrite` failed: Did not find an occurrence of the pattern
  upperRamificationGroup (↥(inertiaGalAdjoin K x)) (stageUniformizer K x) v
in the target expression
  Subgroup.map (inertiaGalAdjoin K x).subtype
      (upperRamificationGroup (↥(inertiaGalAdjoin K x)) (stageUniformizer K x) v) =
    upperRamificationGroup Gal(↥K.carrier⟮x⟯/K.carrier) α v
```

★★**パターンと的が印字上まったく同じ**である。違うのは印字されない
`[Fintype ↥(inertiaGalAdjoin K x)]` インスタンス引数だけで、
`upperRamificationGroup` の型がそれを取っている。

**原因**: 証明の冒頭に
`haveI : Fintype ↥(inertiaGalAdjoin K x) := Fintype.ofFinite _` と書いたため、
的の中では木の canonical instance（`RamificationFiltrationBuild.lean` の
`fintypeInertiaGalAdjoin`）が使われ、`rw` する補題の側では `Fintype.ofFinite _` が
使われて、**構文が一致しない**。

**直し方**: ★**その `haveI` を消す。** 木に instance が既に在るなら書かない。
★一般則: `rw` が「同じに見える項」で失敗したら、**インスタンス引数を疑う**
（`set_option pp.explicit true` で 1 回だけ見ると 5 秒で分かる）。
★これは #126（型に現れるインスタンスは `haveI` では間に合わない）の**裏側**で、
あちらは「足りない」、こちらは「**余計**」である。

## #246 `IsGalois` から `Normal` を `.toNormal` で取り出せない（2026-09-08、同上）

```
error: Invalid field `toNormal`: The environment does not contain `IsGalois.toNormal`, so it is not possible to project the field `toNormal` from an expression
  (InfiniteGalois.normal_iff_isGalois K.carrier⟮(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt⟯).mp hnorm
of type
  IsGalois K.carrier ↥K.carrier⟮(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt⟯
```

**直し方**: フィールド射影ではなく**インスタンス探索**に載せる:

```lean
theorem normal_… : Normal K.carrier ↥L := by
  haveI := isGalois_… ;  infer_instance
```

★`IsGalois` は `Normal` と `Algebra.IsSeparable` を **instance として持つ**が、
`toNormal` という名前の射影は無い。★同じ形は `IsGalois.toIsSeparable` でも起きる。

## #247 `Subgroup.map_comap_eq` は `f.range` を**左**に置く（2026-09-08、同上）

```
error: Tactic `rewrite` failed: Did not find an occurrence of the pattern
  ?a ⊓ ⊤
in the target expression
  ⊤ ⊓ upperRamificationGroup Gal(↥K.carrier⟮x⟯/K.carrier) α v = upperRamificationGroup Gal(↥K.carrier⟮x⟯/K.carrier) α v
```

`rw [Subgroup.map_comap_eq, Subgroup.range_subtype, htop]` のあとに
`inf_top_eq` を書いて落ちた。★`Subgroup.map_comap_eq f H = f.range ⊓ H` なので
正しくは **`top_inf_eq`** である。
★1 秒で直るが、`⊓` の向きは `Submodule` / `Ideal` 版と揃っていないので毎回測ること。

## #248 ★★★「段データが選択に依らない」は `compat` を `M = N` に当てるだけで出る（2026-09-08、同上）

配管ではなく**設計**の記録。`StageFiltration` のような
「`N ≤ M` のとき `S N v · M = S M v`」型の両立性を持つ構成では、
★**`M := N` を代入すると、`N ≤ S N v` と合わせて `S N v = S' N v`（別の選択で作った側）が出る**。

```lean
theorem coe_mul_coe_eq_self {Γ : Type*} [Group Γ] {S N : Subgroup Γ} (h : N ≤ S) :
    ((S : Set Γ) * (N : Set Γ)) = (S : Set Γ) := …   -- 2 行、[propext, Quot.sound]

theorem stage_eq_stage (g g' : StageGenerator K N) (v : ℝ) : g.stage v = g'.stage v :=
  eq_of_coe_mul_coe_eq_coe (le_stage g v) (stage_mul_coe_eq h K g g' le_rfl v)
```

★実測: `Found/PGC/RamificationFiltrationBuild.lean` は「生成元 `x` の選択に依らないことは
**証明していない**」（逸脱 2）と書き、`UnramifiedBaseChangeInvariance.lean` も
「★証明していないこと」に挙げていた。★**その `compat` 自身が独立性を含んでいた。**
★教訓: **「両立性を証明したのに独立性が未証明」と書いてあったら、まず `M = N` を代入する。**

## #250 `Set.mem_image` を `obtain` で開くと、出てくる等式は**β 簡約されていない**ので `rw` が当たらない（2026-09-08、pGC「p 進対数が単数を整数環に写す」）

`rw [← image_… ] at h` のあと `obtain ⟨u, hu, hux⟩ := h` とすると、`hux` の左辺は
**関数適用のまま**（`(fun u => ↑↑u - 1) u`）で残る。ゴールの側は β 簡約済みなので
`rw [hux]` が当たらない。逐語のエラー文:

```
error: Tactic `rewrite` failed: Did not find an occurrence of the pattern
  (fun u => ↑↑u - 1) u
in the target expression
  (↑p ^ r)⁻¹ * padicLog K (↑↑u - 1) = ↑(Multiplicative.toAdd y)
…
hux : (fun u => ↑↑u - 1) u = x
⊢ (↑p ^ r)⁻¹ * padicLog K (↑↑u - 1) = ↑(Multiplicative.toAdd y)
```

★**仮定とゴールが画面上ほぼ同じ字面なのに当たらない**ので、見つけにくい。
直し方は**型を書いて写し直す**だけ（`show`/`simp only []`/`beta_reduce` より短い）:

```lean
have hux' : (((u : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1 = x := hux
rw [hux']
```

★同じことは `Set.BijOn.surjOn` / `Set.image_eq` を経由した `obtain` すべてで起きる。

## #251 `c • (c⁻¹ * y) = y` に `field_simp` を撃つと止まる —— `mul_inv_cancel_left₀` を直に使う（2026-09-08、同上）

`Set.smul_set` の逆像を作るところで出る形。`simp only [smul_eq_mul]` のあとに
`field_simp` を置くと、体でない（`NormedDivisionRing`）ため次で止まる:

```
error: `field_simp` made no progress on the goal
```

`exact mul_inv_cancel_left₀ hc y`（`hc : c ≠ 0`）で 1 行。
★抽象核 `smul_setOf_norm_le : c • {x | ‖x‖ ≤ s} = {y | ‖y‖ ≤ ‖c‖ * s}`
（`Found/PGC/PadicLogIntegers.lean`）はこれで閉じる。

## #252 `Ideal.ramificationIdx` の同定に Dedekind 環の因子分解は要らない —— 第 2 条件はノルム 1 本（2026-09-08、同上）

`Ideal.ramificationIdx_spec (hle : map f p ≤ P ^ n) (hgt : ¬ map f p ≤ P ^ (n+1)) : ramificationIdx f p P = n`
の `hgt` は、離散付値環では **`‖π‖ < 1` だけ**で出る（`IsDedekindDomain.ramificationIdx_eq_*` を通さない）:

```lean
intro hcon
have hmem := hcon (Ideal.mem_span_singleton_self _)   -- p ∈ 𝔪^(e+1)
rw [Ideal.span_singleton_pow] at hmem
have hn := norm_le_of_mem_span_pow K (e + 1) _ hmem   -- ‖p‖ ≤ ‖π‖^(e+1)
rw [h3, hnorm, pow_succ] at hn                        -- ‖π‖^e ≤ ‖π‖^e * ‖π‖
nlinarith [pow_pos hπpos e, norm_pi_lt_one K hπmax]
```

★`map (algebraMap ℤ_[p] 𝒪_K) (maximalIdeal ℤ_[p]) = span {(p : 𝒪_K)}` は
`rw [PadicInt.maximalIdeal_eq_span_p, Ideal.map_span]; simp` で出る（`simp` が像の単集合を潰す）。
★★**mathlib の名前は `IsDiscreteValuationRing.ideal_eq_span_pow_uniformizer` ではなく
`…ideal_eq_span_pow_irreducible`**（`Unknown constant` を撃った）。
`x = u * ϖ^n` が欲しいだけなら `IsDiscreteValuationRing.eq_unit_mul_pow_irreducible` の方が直接。

## #253 `σ • x` と `σ x` は defeq だが **`rw` は当たらない** —— `exact`/`show` に替える（2026-09-08、pGC `HasCoherentFunctional`）

`MulSemiringAction` 由来の `•` は適用と defeq だが、**`rw` は構文照合**なので落ちる。
`coe_levelPi_apply (g) (b) : ↑((levelPi K L g) b) = g ↑b` を `↑(levelPi K L' γ • a)` に当てようとすると:

```
error: Tactic `rewrite` failed: Did not find an occurrence of the pattern
  ↑(((levelPi ?K ?L) ?g) ?b)
in the target expression
  ↑((levelPi K L') γ • a) = ↑a
```

同じ罠は `exact` の**向きを間違えたとき**に別の顔で出る（★こちらは defeq が効いているので
`.symm` を外すだけで通る。「型が合わない」ように見えて実は向きの問題）:

```
error: Type mismatch
  Eq.symm (coe_levelPi_apply K L' γ (avgSum (levelRes K L L' h).ker x))
has type
  γ ↑(avgSum (levelRes K L L' h).ker x) = ↑(((levelPi K L') γ) (avgSum (levelRes K L L' h).ker x))
but is expected to have type
  ↑((levelPi K L') γ • avgSum (levelRes K L L' h).ker x) = γ ↑(avgSum (levelRes K L L' h).ker x)
```

★直し方は 2 つ。(i) `rw` をやめて `exact <補題>`（defeq で通る）。
(ii) 目標側を `show` で**適用の形に書き換えてから** `rw` する:
`show ((levelPi K L' γ (levelIncl K L L' h b) : L') : K.closure) = _`。
★`rw [← coe_levelPi_apply K L' γ a]` と**逆向きに使う**手もある（`γ ↑a` を `↑(… a)` に戻す）。

## #254 データを `have` で置くと本体を忘れる —— `show` が「defeq でない」と言い出す（2026-09-08、同上）

`have e : N ≃ N := { toFun := …, invFun := …, … }` としてから `Fintype.sum_bijective e e.bijective`
に渡すと、最後の点ごとの等式で:

```
error: 'show' tactic failed, pattern
  g • ↑n • a = (g * ↑n * g⁻¹) • g • a
is not definitionally equal to target
  g • ↑n • a = ↑(e n) • g • a
```

`have` は **Prop でなくても本体を捨てる**（`e` は不透明な局所仮定になる）。
★直し方: 関数を `refine` の中に**直に書く**（`refine Fintype.sum_bijective (fun n : N => ⟨g * n * g⁻¹, …⟩) ?_ _ _ ?_`）。
全単射性は `Function.bijective_iff_has_inverse` に逆写像を渡すのが軽い（`Equiv` を組み立てなくてよい）。

## #255 `Finset.sum_congr rfl (fun n _ => h n n.2)` の `n.2` が **G の元に潰れる**（2026-09-08、同上）

`∑ n : ↥N, (n : G) • a` の項に `Finset.sum_congr` を当てると、ラムダの `n` が
**強制で `G` に coerce されてから**射影を取ろうとして落ちる:

```
error: Invalid projection: Projection operates on types of the form `C ...` where C is a constant. The expression
  n
has type `G` which does not have the necessary form.
```

★直し方: 先に `have hn : ∀ n : N, (n : G) • a = a := fun n => h (n : G) n.2` と**外へ出して**から
`rw [Finset.sum_congr rfl (fun n _ => hn n)]`。

## #256 `F⟮α⟯` 記法は**スコープ外だと 2 つのエラーに割れる**（2026-09-08、同上）

`x ∈ K.carrier⟮α⟯` と書くと、記法が入っていない環境では `∈` までで切れて:

```
error(lean.synthInstanceFailed): failed to synthesize instance of type class
  Membership K.closure Type
error: expected token
```

★1 つ目のエラーが「`Membership _ Type`」——**右辺が `Type` になっている**のが合図である
（`K.carrier` そのものを集合と読んでいる）。
★直し方: `IntermediateField.adjoin K.carrier {α}` と書く（記法に依存しない）。
`IntermediateField.adjoin_simple_le_iff : F⟮α⟯ ≤ K ↔ α ∈ K` はそのまま当たる（定義が同じ）。

## #257 `.lean`/`.md` を `cat > f <<'EOF'` で書くと**この環境では壊れる**（2026-09-08、同上）

Bash の PreToolUse フックがコマンドを書き換えるため、引用符付き heredoc でも:

```
/usr/bin/bash: -c: line 188: unexpected EOF while looking for matching `''
```

（★ファイルは**作られない**ので実害は往復 1 回ぶんだが、`>` が先に効く形だと切り詰めが起きうる。）
★直し方: **Write/Edit ツールで書く**。★`sed -i` のような 1 行編集は通る。

## #258 並行セッションの `lake build` は **`no such file or directory` で落ちる** —— 自分のエラーではない（2026-09-08、pGC `HasCoherentFunctional`）

同じワークツリーで別の agent が `lake build` していると、olean の書き込みが衝突して落ちる:

```
✖ [3199/3232] Building ABC3.Found.PGC.LubinTateReciprocityLimitCompat (8.5s)
error: no such file or directory (error code: 4294963238)
  file: D:\Math_ABC3\lean\.lake\build\lib\lean\ABC3\Found\PGC\LubinTateReciprocityLimitCompat.olean
```

`leanfile.mjs` でも同じ原因で別の顔になる（★依存の olean が**消えている**瞬間に読む）:

```
ABC3/Found/PGC/CoherentFunctional.lean:1:0: error: object file
'…\ABC3\Found\PGC\LubinTateGeneralUniqueness.olean' of module
ABC3.Found.PGC.LubinTateGeneralUniqueness does not exist
```

★**合図**: 落ちている宣言が**自分の持ち場と無関係**で、しかも
`build.mjs` の要約が `error 0 / sorry 0` のまま「★失敗」になる。
★直し方: **もう一度同じコマンドを打つだけ**（2026-09-08 に 2 回とも 1 回の再実行で通った。
198.6 秒 → 116.4 秒、379.3 秒 → 75.3 秒）。★ファイルを直さないこと。

## #259 `pow_lt_pow_left` は改名された + `zero_le` は NNReal では**関数ではない**（2026-09-08、同上）

```
error(lean.unknownIdentifier): Unknown identifier `pow_lt_pow_left`
```

★正しい名前は **`pow_lt_pow_left₀ (hab : a < b) (ha : 0 ≤ a) (hn : n ≠ 0) : a ^ n < b ^ n`**
（`Mathlib/Algebra/Order/GroupWithZero/Basic.lean`）。順序付きモノイド版は `pow_lt_pow_left'`。
★続けて `(zero_le _)` と書くと NNReal では次で落ちる（`zero_le : 0 ≤ a` は**引数を取らない**）:

```
error: Function expected at
  zero_le
but this term has type
  0 ≤ ?m.968
```

★直し方: `pow_lt_pow_left₀ hδ2 zero_le hd.ne'`。

## #260 宣言を `Found/` から `Skeleton/` へ移すと **`check.mjs` G1 が `.src` を新たに要求する**（2026-09-08、pGC `FilteredGroup` の移設）

`Found/PGC/FilteredGroup.lean` にあった `structure FilteredGroup.Iso` /
`def FilteredGroup.OuterIso` を `Skeleton/PGC/Setup.lean` へ移したところ、
★**中身を 1 文字も変えていないのに** `node tools/check.mjs --ledger --brief` が 2 件増えた:

```
NG  lean\ABC3\Skeleton\PGC\Setup.lean:215
      G1 出典が無い: `FilteredGroup.Iso.src : ABC3.Meta.Source` を書く
NG  lean\ABC3\Skeleton\PGC\Setup.lean:234
      G1 出典が無い: `FilteredGroup.OuterIso.src : ABC3.Meta.Source` を書く
```

★G1 は `Skeleton/` と `Interface/` の宣言にだけ `.src` を要求する（`Found/` には要求しない）。
★**移設は「移すだけ」では終わらない**——移す前に、移す本の bucket が変わるかを見て、
変わるなら `.src` を同時に書くこと。

★原典が定義を読者に委ねている（我々自身の定式化である）場合は、bare な `"Definition 2.3"` に
せず `item := "Definition 2.3 (FilteredGroup.Iso)"` と注記を付ける
——bare だと G9（非空虚性の対照）の対象になる。

★逆向きの注意: `Interface/` へ移すと今度は **G2**（`check.mjs:1041`）が
`X.nonvacuous` か `X.waiting` を要求する。★`structure` を移す先を決める前に、
G1 / G2 / G8 / G9 のどれが新たに掛かるかを数えること。

## #261 node で `.lean`/`.md` を `latin1` で読むと、UTF-8 の文字列リテラルと**照合が必ず外れる**（2026-09-08、同上）

CRLF を保ったまま行を消す・入れ替えるスクリプトを書くとき、改行を数える都合で
`fs.readFileSync(p, 'latin1')` としたくなる。★しかしその文字列に対して
スクリプト中の日本語リテラル（node は UTF-8 で読む）を `includes` すると**必ず false** になる:

```
NG: docstring のアンカーが無い
```

★ASCII だけの検算（`/-! ## Corollary 3.1 -/` など）は通ってしまうので、
**「一部だけ通って一部だけ落ちる」**という分かりにくい形で出る。
★直し方: **読むのも書くのも `utf8`**。`\r` は普通の文字なので CRLF は utf8 でもそのまま保たれる
（`(s.match(/\r\n/g)||[]).length` で前後を数えて確かめる）。
★書き出しは `fs.writeFileSync(p+'.tmp', Buffer.from(out,'utf8'))` → `fs.renameSync` の順にする
（失敗しても元ファイルが無傷）。

## #262 **型同義語**（`def T := ℚ_[p]`）を実引数にすると、暗黙引数の単一化が `whnf` で止まる（2026-09-08、pGC Proposition 2.2）

`Check/PGC/Prop12Degenerate.lean` の `TwistedQp p : Type := ℚ_[p]`（体構造だけ捻った型同義語）を
`{K K' : PAdicLocalField p}` が暗黙の補題へ渡すと、**証明項を書いただけで**次が出る:

```
error: (deterministic) timeout at `whnf`, maximum number of heartbeats (200000) has been reached
```

★エラーは補題の適用行ではなく **`theorem` の宣言行（`7:0`）** に出るので、
「文が重いのか証明が重いのか」が分からない。★実測では**文だけなら通る**（`sorry` を置くと 1 秒）。

★原因: `?K.carrier =?= TwistedQp p` の右辺が `ℚ_[p]` に簡約できてしまうため、
単一化器が `selfField p` 側の `ℚ_[p]` と往復して探索が爆発する。

★直し方は **`maxHeartbeats` を上げることではない**（上げても通るが遅い）。
**暗黙引数を名前で固定する**:

```
NG:  intKbar_transport_galContinuousMulEquiv (twistedAlgEquiv p)
OK:  intKbar_transport_galContinuousMulEquiv (K := twistedField p) (K' := selfField p)
       (twistedAlgEquiv p)
```

★実測: `set_option maxHeartbeats` **なしで通る**（12 秒 → 12 秒、`leanfile.mjs`）。
★同じ型同義語を**文の側**でも使うなら、そこにも `(K := …) (K' := …)` を書くこと
（`galContinuousMulEquiv (K := twistedField p) (K' := selfField p) β`）。
★別の顔で出ることもある —— 実引数の型が確定していないと先に

```
error: Application type mismatch: The argument
  twistedAlgEquiv p
has type
  TwistedQp p ≃ₐ[ℚ_[p]] ℚ_[p]
but is expected to have type
  (twistedField p).carrier ≃ₐ[ℚ_[p]] PAdicLocalField.carrier (?m.43 φ g x)
```

が出る。★これも同じ直し方（名前で固定）で消える。

## #263 `Module R (UniformSpace.Completion α)` が付かないのは **`UniformContinuousConstSMul` が無い**から（2026-09-08、pGC §3 の `ℂ_K`）

`CompKbar K = closureCompletion K = UniformSpace.Completion K.closure` に対し、

```
error(lean.synthInstanceFailed): failed to synthesize instance of type class
  Module ℚ_[p] (closureCompletion K)
```

★`Algebra ℚ_[p] K.closure` も `Module ℚ_[p] K.closure` も**在る**（`inferInstance` が通る）。
足りないのは 1 つだけ:

```
error(lean.synthInstanceFailed): failed to synthesize instance of type class
  UniformContinuousConstSMul ℚ_[p] K.closure
```

`UniformSpace.Completion.instModule`（`Topology/Algebra/GroupCompletion.lean:204`）が
`[UniformContinuousConstSMul R α]` を要求しているためである
（`UniformSpace.Completion.algebra`（`Topology/Algebra/UniformRing.lean:198`）も同じ）。

★★**直し方（菱形を作らない側）**: 既存の `Algebra` を**再利用**して `NormedAlgebra` を足す。
`IsBoundedSMul.toUniformContinuousConstSMul`（`Topology/MetricSpace/Algebra.lean:157`、priority 100）が
残りを繋ぐので、これ 1 本で `Module` / `Algebra` / `IsScalarTower` が全部付く:

```
noncomputable scoped instance closureNormedAlgebraQp (K : PAdicLocalField p) :
    NormedAlgebra ℚ_[p] K.closure :=
  { (inferInstance : Algebra ℚ_[p] K.closure) with
    norm_smul_le := fun r x => by
      rw [Algebra.smul_def, norm_mul, IsScalarTower.algebraMap_apply ℚ_[p] K.carrier K.closure,
        norm_algebraMap_closure, ABC3.Found.PGC.norm_algebraMap] }
```

★★**やってはいけない側**: `((algebraMap S A).comp (algebraMap R S)).toAlgebra` で
`Algebra R (Completion α)` を**新しく**立てること。`SMul` が
`UniformSpace.Completion.instSMul` と割れて、症状は `show` の失敗として出る:

```
error: 'show' tactic failed, pattern
  σ •
      @HSMul.hSMul ℚ_[p] (CompKbar K) (CompKbar K)
        (@instHSMul ℚ_[p] (CompKbar K) (UniformSpace.Completion.instSMul ℚ_[p] K.closure)) c z =
    c • σ • z
is not definitionally equal to target
  σ •
      @HSMul.hSMul ℚ_[p] (CompKbar K) (CompKbar K) (@instHSMul ℚ_[p] (CompKbar K) DistribMulAction.toDistribSMul.toSMul)
        c z =
    (RingHom.id ℚ_[p]) c • σ • z
```

★同じ回で出た小物 2 つ:

```
error: Type mismatch
  IsScalarTower.of_algebraMap_eq fun x ↦ rfl
has type
  IsScalarTower ?m.10 ?m.10 ?m.22
but is expected to have type
  IsScalarTower ℚ_[p] K.carrier (CompKbar K)
```

→ `fun _ => rfl` は何の情報も与えないので `(R := …) (S := …) (A := …)` を書く。
★ただし上の `NormedAlgebra` を足した後は `inferInstance` で付くので、そもそも要らなくなる。

```
error(lean.unknownIdentifier): Unknown constant `MonoidHom.fst_apply`
```

→ そんな simp 補題は無い。`(MonoidHom.fst M N) (a, b) = a` は `rfl` なので、
`simp only [...]` の後に `rfl` を 1 行置く。

## #264 `Finset.sum_erase _ (h : f a = 0)` は **`f` と `a` を取り違えて単一化する**（2026-09-08、pGC Ax–Sen–Tate の抽象核）

`∑ i ∈ s.erase i0, d i • w i = ∑ i ∈ s, d i • w i` を出したくて
`Finset.sum_erase`（`(s : Finset α) (h : f a = 0) : ∑ x ∈ s.erase a, f x = ∑ x ∈ s, f x`）に
型注釈だけ付けて渡すと落ちる:

```
rw [Finset.sum_erase _ (by rw [hdi0, zero_smul] : d i0 • w i0 = 0)]
```

```
error: Tactic `rewrite` failed: Did not find an occurrence of the pattern
  ∑ x ∈ Finset.erase ?m.545 (w i0), d i0 • x
in the target expression
  ∑ i ∈ s.erase i0, d i • w i = 0
```

★エラー文の `∑ x ∈ Finset.erase ?m (w i0), d i0 • x` が症状を全部語っている ——
`h` の型 `d i0 • w i0 = 0` から `f := fun x => d i0 • x`、`a := w i0` と読まれた。
`f a = 0` は `f` と `a` の**両方が高階の未知数**なので、単一化は左から一番浅い分解を選ぶ。

★★**直し方**: `f` を名前付き引数で釘付けにする（`a` は `f` が決まれば決まる）。

```
rw [Finset.sum_erase (f := fun i => d i • w i) s (by rw [hdi0, zero_smul])]
```

★同族: `Finset.prod_erase` / `Finset.sum_subset` / `Finset.sum_congr` の
「`h : f a = 0` から `f` を推論させる」形は全部これに当たる。

## #265 `TensorProduct.smul_tmul'` は **名前から想像する向きの逆**（2026-09-08、同上）

```
theorem TensorProduct.smul_tmul' (r : R') (m : M) (n : N) : r • m ⊗ₜ[R] n = (r • m) ⊗ₜ n
```

★**左辺が `r • (m ⊗ₜ n)`**（スカラーが外）で、**右辺が `(r • m) ⊗ₜ n`**（スカラーが中）である。
`'` が付いているので「`smul_tmul` の逆向き」と思って `← smul_tmul'` を書くと落ちる:

```
error: Tactic `rewrite` failed: Did not find an occurrence of the pattern
  (?r • ?m) ⊗ₜ[?m.94] ?n
in the target expression
  (((compRep K).tprod ρ) σ) (a • x ⊗ₜ[ℚ_[p]] v) = (σ • a) • (((compRep K).tprod ρ) σ) (x ⊗ₜ[ℚ_[p]] v)
```

★★**紛らわしさの正体**: `a • x ⊗ₜ[ℚ_[p]] v` という印字は
`a • (x ⊗ₜ v)` と `(a • x) ⊗ₜ v` の**どちらにも読める**（`•` の方が `⊗ₜ` より強く結合するので
実際は前者）。だから「もう在る形」を目で判定できない。★エラー文の `pattern` の側を見ること
——`←` を付けたときに表示される pattern は**補題の右辺**である。

★**直し方**: `a • (x ⊗ₜ v)` を潰したいなら `←` を**付けない**。
`simp only [TensorProduct.smul_tmul', …, smul_eq_mul, smul_mul']` の 1 行で、
半線型性 `σ • (a • z) = (σ • a) • (σ • z)` の `tmul` の場合が閉じる。

## #266 `Cardinal.nat_lt_aleph0` の後継 `Cardinal.natCast_lt_aleph0` は **引数が暗黙**（2026-09-08、同上）

```
warning: `Cardinal.nat_lt_aleph0` has been deprecated: Use `Cardinal.natCast_lt_aleph0` instead
```

に従って機械的に置き換えると次で落ちる:

```
error: Function expected at
  Cardinal.natCast_lt_aleph0
but this term has type
  ↑?m.88 < Cardinal.aleph0
```

→ `Cardinal.nat_lt_aleph0 n` は明示引数だったが、後継は暗黙引数なので
**`_` ごと消す**（`Cardinal.natCast_lt_aleph0` と書く）。
★deprecation の指示は名前しか直してくれない。★引数の explicit/implicit まで変わる例である。


## #267 `rw` に渡す `have` は**ゴールと同じ coercion の書き方**で述べる —— `MulEquiv.toMonoidHom e` は `↑e` に当たらない（2026-09-08、濾過つき `IntKbarRecoverable`）

`structure FilteredGroup.Iso` のフィールド
`map_Gv : ∀ v, Subgroup.map ↑equiv.toMulEquiv (A.Gv v) = B.Gv v` に対し、
`Subgroup.map_map` のあと「合成が恒等」を `have` で用意して `rw` すると落ちる:

```
error: Tactic `rewrite` failed: Did not find an occurrence of the pattern
  f.equiv.symm.toMonoidHom.comp f.equiv.toMonoidHom
in the target expression
  Subgroup.map ((↑f.equiv.symm.toMulEquiv).comp ↑f.equiv.toMulEquiv) (A.Gv v) = A.Gv v
```

`have hid : (MulEquiv.toMonoidHom e.symm).comp (MulEquiv.toMonoidHom e) = MonoidHom.id _` と
書くと `e.toMonoidHom` に elaborate されるが、ゴールに在るのは `↑e`（`MonoidHomClass`
経由の coe）である。**両者は defeq だが構文が違う**ので `rw` が当たらない。
**直し方**: `have` の型を **`((e : A ≃* B) : A →* B)` の形（`↑` と同じ書き方）**で述べる。

```lean
have hid : ((f.equiv.symm.toMulEquiv : B.G ≃* A.G) : B.G →* A.G).comp
    ((f.equiv.toMulEquiv : A.G ≃* B.G) : A.G →* B.G) = MonoidHom.id A.G :=
  MonoidHom.ext (fun x => f.equiv.symm_apply_apply x)
rw [← f.map_Gv v, Subgroup.map_map, hid, Subgroup.map_id]
```

★これで一発で通った（`node tools/leanfile.mjs`、第 1073）。
★#96「`set` で束縛した局所定義は `rw` のパターンに当たらない」と同じ
「**欲しい形で型を書き直す**」の coercion 版である。

## #268 構造体を返す `def` の射影は **instance 探索を通さない**（型検査は defeq で通るのに）—— 型注釈つきの薄い `def` を 1 本挟む（2026-09-08、濾過つき `IntKbarRecoverable`）

`noncomputable def pgcFilteredGroup (K : PAdicLocalField p) : FilteredGroup :=
filteredGroupOf (ramificationFiltration p) K` と置き、
`α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K')` から
`(α.equiv.toMulEquiv g) • x` と書くと落ちる:

```
error(lean.synthInstanceFailed): failed to synthesize instance of type class
  HSMul (pgcFilteredGroup K').G (Obj K') ?m.41
```

`(pgcFilteredGroup K').G` は `K'.absGal` と defeq である——同じ木で
`example (α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K')) :
ContinuousMulEquiv K.absGal K'.absGal := α.equiv` は**通る**。
落ちるのは **instance 探索が `def` を delta 展開しない**からで、型検査の失敗ではない。
**直し方**: 型注釈つきの薄い `def` を 1 本挟み、以後はそれだけを使う。

```lean
def filteredIsoEquiv {K K' : PAdicLocalField p}
    (α : FilteredGroup.Iso (pgcFilteredGroup K) (pgcFilteredGroup K')) :
    ContinuousMulEquiv K.absGal K'.absGal := α.equiv
```

★挟んだあとは `filteredIsoEquiv (filteredIsoRefl _) = ContinuousMulEquiv.refl _` など
4 本が**すべて `rfl`** で通った（1 層の射影だから。2 層は #59 で止まる）。
★既出節「構造体を返す `def` の射影ではインスタンスが見つからない」（`SSCurve.ext`）と
同じ罠だが、そこにはエラー文が 1 つも引用されていなかったので逐語で足した。

## #269 `ring` は**非可換 `[Ring A]`** でも `ring_nf made no progress` になる —— 抽象核を `AddCommGroup` + `nsmul` に落とすと消える（2026-09-08、pGC Ax の補題）

**症状**: 「Multiset の各項から定数を引いた和」の抽象核を `[Ring A]` で書いた:

```lean
theorem multisetSum_map_const_sub {A : Type*} [Ring A] (s : Multiset A) (x : A) :
    (s.map (fun r => x - r)).sum = (Multiset.card s : A) * x - s.sum := by
  induction s using Multiset.induction with
  | cons a t ih => rw [...]; push_cast; ring
```

```
error: `ring_nf` made no progress on the goal
```

**原因**: `ring` は **`CommSemiring`** を要求する。`[Ring A]` には `+` も `*` も在るので
#178（`CommMonoid` に `+` が無い）とは別の理由だが、**エラー文の顔は同じ**である。
`(Multiset.card s : A) * x` の `natCast` が `x` と可換であることを `ring` は使えない。

**直し方（★一般化すると消える）**: 結論を `nsmul` で述べ、`[AddCommGroup A]` まで弱めて `abel`。

```lean
theorem multisetSum_map_const_sub {A : Type*} [AddCommGroup A] (s : Multiset A) (x : A) :
    (s.map (fun r => x - r)).sum = Multiset.card s • x - s.sum := by
  induction s using Multiset.induction with
  | empty => simp
  | cons a t ih =>
      rw [Multiset.map_cons, Multiset.sum_cons, ih, Multiset.sum_cons, Multiset.card_cons,
        succ_nsmul]
      abel
```

★これで `#print axioms` が `[propext, Quot.sound]` だけになった（`Classical.choice` すら要らない）。
★使う側（体）では `nsmul_eq_mul` で `(card : A) * x` に戻せばよい。

## #270 群作用を抽象核に食わせるときは **`M →+ M` に名前を付けて `@[simp] _apply` を 1 本置く** —— 無名の `DistribSMul.toAddMonoidHom` は `rw` が当たらない（2026-09-08、同上）

**症状**: 抽象核が `act : G → M →+ M` を取るので、具体層で
`fun σ => DistribMulAction.toAddMonoidHom (CompKbar K) σ` を直接渡したところ、まず

```
warning: `DistribMulAction.toAddMonoidHom` has been deprecated: Use `DistribSMul.toAddMonoidHom` instead
```

が出て、さらにゴールに現れた項に `smul_coe_closureCompletion`（`σ • ↑x = ↑(σ • x)`）が
**当たらず**、`rfl` で閉じようとすると

```
Tactic `rfl` failed: The left-hand side
  (DistribMulAction.toAddMonoidHom (CompKbar K) σ) ↑x0 - ↑x0
is not definitionally equal to the right-hand side
  ↑(σ • x0) - ↑x0
```

**原因**: `rw` は**構文**で当てるので、`(toAddMonoidHom … σ) w` の中の `σ • w` は見えない。
一方 `↑(σ • x0)`（`K̄` 側の作用）と `σ • ↑x0`（完備化側の作用）は
**定義的には等しくない**（後者は一様連続延長）ので `rfl` も通らない。

**直し方**: 名前付きの薄い `def` と `@[simp]` の `_apply` を置く。`_apply` は `rfl` で通る。

```lean
noncomputable def smulAddHom (K : PAdicLocalField p) (σ : K.absGal) :
    CompKbar K →+ CompKbar K :=
  DistribSMul.toAddMonoidHom (CompKbar K) σ

@[simp] theorem smulAddHom_apply (K : PAdicLocalField p) (σ : K.absGal) (w : CompKbar K) :
    smulAddHom K σ w = σ • w := rfl
```

★`noncomputable` を忘れると

```
error(lean.dependsOnNoncomputable): failed to compile definition, consider marking it as
'noncomputable' because it depends on 'UniformSpace.Completion.instAddMonoid',
which is 'noncomputable'
```

★以後 `rw [smulAddHom_apply]` で `σ • w` が表に出るので、`smul_coe_closureCompletion` が当たる。

## #271 ★訂正: `lean/ABC3/Found.lean` は **2026-09-08 時点で LF** —— #141 / #185 の「CRLF」は古い（2026-09-08、同上）

**症状**: #141 / #185 は「`Found.lean` は CRLF だから Python は `'rb'` で読め」と言う。
実測（raw byte を数えた）:

```
$ node -e "const b=require('fs').readFileSync('D:/Math_ABC3/lean/ABC3/Found.lean');
  let cr=0,lf=0; for(const c of b){if(c===13)cr++;if(c===10)lf++;} console.log('CR',cr,'LF',lf)"
CR 0 LF 1836
```

`lean/ABC3/Found/PGC/AxSenTate.lean` も `CR 0 LF 505`。
★**CRLF で書き戻すと 1 ファイル丸ごと差分になる**（#141 が警告していたのと逆向きの事故）。

**★測り方**: `node tools/eol-audit.mjs --tracked` の一覧に**載っていなければ LF** である。
★**`grep -c $'\r' FILE` を使ってはならない** —— この木の Bash ツールでは `$'\r'` が
展開されず**空パターンになって全行にマッチ**し、LF のファイルが「全行 CRLF」に見える
（実測でそう誤読した）。`cat -A` は #185 のとおり CR を落とす。
★raw byte を数えるか `eol-audit.mjs` を使うこと。

## #272 作用の定義を展開するとき `rw [smul_..._def]` は**最初の 1 箇所しか当たらない** —— `simp only` にする（2026-09-08、pGC Sen の補題）

**症状**: `σ • r` を `σ r` に開いて `map_div₀` を当てようとすると、
ゴールの**左辺だけ**が開いて右辺の `σ • (a / b)` が残り、次の `rw` が落ちる。

```
error: Tactic `rewrite` failed: Did not find an occurrence of the pattern
  ?f (?a / ?b)
in the target expression
  (Polynomial.aeval (σ r)) g / (Polynomial.aeval (σ r)) (Polynomial.derivative (minpoly K.carrier x)) =
    σ • ((Polynomial.aeval r) g / (Polynomial.aeval r) (Polynomial.derivative (minpoly K.carrier x)))
```

**原因**: `rw` は**最左最内の 1 インスタンス**だけを書き換える。
`smul_closure_def : σ • x = σ x`（`Found/PGC/AbsClosureModules.lean:260`）は
`σ • _` の形なので、`w (σ • r) = σ • w r` のような**両辺に `•` が出る同変性の証明**では
片側しか開かない。

**直し方**: 定義展開は `simp only` にする。

```lean
×  rw [smul_closure_def, map_div₀, ...]
○  simp only [smul_closure_def]
   rw [map_div₀, Polynomial.aeval_algHom_apply σ r g, ...]
```

★同じ罠は `smul_eq_mul` / `Algebra.smul_def` / `Module.End` 系の `_def` すべてに起きる。
★**「`_def` は `simp only`、補題は `rw`」**と覚えるとよい。

## #273 `Polynomial.derivative_map` の向き —— `←` を付けると当たらない（2026-09-08、同上）

**症状**:

```
error: Tactic `rewrite` failed: Did not find an occurrence of the pattern
  Polynomial.map ?f (Polynomial.derivative ?p)
in the target expression
  (Polynomial.aeval r) (Polynomial.derivative f) =
    Polynomial.eval r (Polynomial.derivative (Polynomial.map (algebraMap K.carrier K.closure) f))
```

**原因**: `Polynomial.derivative_map : (p.map f).derivative = p.derivative.map f` は
**`derivative (map …)` → `map (derivative …)`** の向きである。
`← derivative_map` は逆向きのパターン `map f (derivative p)` を探すので、
`derivative (map f p)` しか無いゴールでは当たらない。

**直し方**: `rw [Polynomial.derivative_map, Polynomial.eval_map, ← Polynomial.aeval_def]`
（`aeval r q = eval r (q.map (algebraMap ..))` へ落とす定型）。

## #274 `Lagrange.coeff_eq_sum` は **Euler の等式 `Σ_r g(r)/f′(r) = 1` の最短路**（2026-09-08、同上）

**在庫**: `Lagrange.coeff_eq_sum`（`LinearAlgebra/Lagrange.lean:495`）

```
theorem coeff_eq_sum (hvs : Set.InjOn v s) {P : Polynomial F} (hP : P.degree < #s) :
    P.coeff (#s - 1) = ∑ i ∈ s, (P.eval (v i)) / ∏ j ∈ s.erase i, (v i - v j)
```

★節点 `v = id` を `minpoly` の根の集合に取り、`Π_{j≠r}(r−j) = f′(r)`
（`Polynomial.eval_multiset_prod_X_sub_C_derivative`, `Algebra/Polynomial/Derivative.lean:680`）
と組むと、**モニックで次数 `n−1` の `g` について `Σ_r g(r)/f′(r) = 1`** が出る。

★★**`Module.Basis.traceDual_powerBasis_eq` 経由より軽い** —— そちらは `PowerBasis K L`、
つまり**中間体 `K(x)` を作る**ことを要求し、`lean-idioms.md` #59/#69 の境界に当たる。
`Lagrange.coeff_eq_sum` なら**根の重複集合の上だけ**で閉じる。

★`Multiset` ↔ `Finset` の往復は分離性（標数 0 ⇒ `PerfectField.ofCharZero` ⇒
`Algebra.IsSeparable.isSeparable` ⇒ `Polynomial.nodup_roots`）で
`Multiset.dedup_eq_self.mpr` / `Multiset.toFinset_val` / `Finset.erase_val` を使う。
`(s.map h).sum = ∑ r ∈ s.toFinset, h r` は `rw [← htval]; rfl` で閉じる。

## #275 `Subgroup.map` に `MulEquiv` をそのまま渡すと落ちる —— 型上昇は `(Φ : Γ →* Γ')` と**書く**（2026-09-08、pGC 分岐濾過の自然性）

**症状**: 抽象核を `Φ Ψ : Γ ≃* Γ'` で書き、結論に `Subgroup.map Φ H = Subgroup.map Ψ H` と書いた。

```
error: Application type mismatch: The argument
  Ψ
has type
  Γ ≃* Γ'
but is expected to have type
  Γ →* ?m.24
in the application
  Subgroup.map Ψ
```

**原因**: `Subgroup.map` は `MonoidHom` を取る。`MulEquiv → MonoidHom` の coe は
**期待型が決まっているときだけ**入る。`Subgroup.map Φ H = RHS` のように RHS が
先に型を決めてくれる書き方（`Interface` の `IsNaturalFiltration` はこの形）では通るが、
`def` の中で左から順に読まれると `?m` が残って落ちる。

**直し方**: `(Φ : Γ →* Γ')` と**型を書く**。`_ →* _` は解けないことがある（#277）。

```lean
theorem map_eq_map_of_conj (Φ Ψ : Γ ≃* Γ') (H : Subgroup Γ) [H.Normal]
    (h : ∀ g : Γ, ∃ c : Γ', Φ g = c * Ψ g * c⁻¹) :
    Subgroup.map (Φ : Γ →* Γ') H = Subgroup.map (Ψ : Γ →* Γ') H := …
```

★**副作用**: こう書くと `↑Φ x` と `Φ x` が別の項として印字されるので、
`rw [hc]`（`hc : Φ x = …`）が当たらなくなる:

```
error: Tactic `rewrite` failed: Did not find an occurrence of the pattern
  Φ x
in the target expression
  ↑Φ x ∈ Subgroup.map (↑Ψ) H
```

★直前に `show Φ x ∈ Subgroup.map (Ψ : Γ →* Γ') H` を 1 行置くと通る（defeq なので `show` で降りる）。

★**さらに**: `Subgroup.mem_map_equiv` はこの形に**当たらない**:

```
error: Tactic `rewrite` failed: Did not find an occurrence of the pattern
  ?m.83 ∈ Subgroup.map (MulEquiv.toMonoidHom ?m.81) ?m.82
in the target expression
  y ∈ Subgroup.map (↑Φ) (⨅ N, A ↑N) ↔ y ∈ ⨅ N', A' ↑N'
```

⇒ `Subgroup.mem_map_equiv` は諦めて `le_antisymm` ＋ `rintro _ ⟨x, hx, rfl⟩` で書く方が速い。

## #276 `rintro _ ⟨x, hx, rfl⟩` で降りたゴールは **`Set` の膜が 1 枚残る** —— `simp only [SetLike.mem_coe, …]` を先頭に置く（2026-09-08、同上）

**症状**: `Subgroup.map f X ≤ Y` を `apply le_antisymm` → `rintro _ ⟨x, hx, rfl⟩` で降ろし、
ゴールに `Subgroup.mem_iInf` を当てようとした。

```
error: `simp` made no progress
```

**原因**: `≤` は `SetLike` 経由で `↑A ⊆ ↑B` に展開されるので、
降りたゴールは `x ∈ ↑(⨅ …)`（**`Set` のメンバーシップ**）であって
`x ∈ (⨅ …)`（`Subgroup` のメンバーシップ）ではない。`Subgroup.mem_iInf` は後者にしか当たらない。

**直し方**: `simp only [SetLike.mem_coe, Subgroup.mem_iInf, Subtype.forall]`。
★**仮定側**（`… at hx`）では `SetLike.mem_coe` が不要なことがあり、
linter が `This simp argument is unused: SetLike.mem_coe` と言う。**ゴール側だけ**に付ける。

## #277 `(e : _ ≃* _) : _ →* _` の二重アスクリプションは**引数位置では解けない**（2026-09-08、同上）

**症状**: `refine map_upperRamificationGroup_eq (f := ((e : _ ≃* _) : _ →* _)) hsurj ?_ v`

```
error: Type mismatch
  stageInertiaEquiv α x
has type
  ↥(inertiaGalAdjoin K x) ≃* ↥(inertiaGalAdjoin K' ((extendToClosure α) x))
but is expected to have type
  ?m.268 →* ?m.269
```

さらに `hsurj` の側でも:

```
error: Type mismatch
  MulEquiv.surjective (stageInertiaEquiv α x)
has type
  Function.Surjective ⇑(stageInertiaEquiv α x)
but is expected to have type
  Function.Surjective ⇑?m.290
```

**直し方（2 つ、どちらも実測で効いた）**:

1. ★**`f` を渡さない**。`refine map_upperRamificationGroup_eq ?_ ?_ v` にすると
   結論とゴールの単一化で `f` が決まり、`· exact e.surjective` が defeq で通る。
2. ★どうしても書くなら **型を全部書く**（`_` を 1 つも使わない）。
   `((e.symm : A ≃* B) : A →* B)` の `A`・`B` を明示すると通った。
   ★`.comp` の中（`f.comp g = h.comp k` の形）では周りが型を決めるので `_ ≃* _` でも解ける。

## #278 `letI := (foo …).normal` はインスタンス探索に**引っかからない** —— 型を書く（2026-09-08、同上）

**症状**: 構造体フィールドから instance を取り出して `letI` で置いた:

```lean
letI := (stageGeneratorMap α ᾱ hfwd hN g).normal
```

```
error(lean.synthInstanceFailed): failed to synthesize instance of type class
  Normal K'.carrier ↥K'.carrier⟮ᾱ g.gen⟯
```

**原因**: フィールドの型は `Normal K'.carrier ↥K'.carrier⟮(stageGeneratorMap …).gen⟯` であり、
`.gen = ᾱ g.gen` は `rfl` だが **`stageGeneratorMap` は reducible ではない**ので、
インスタンス探索（reducible 透明度）は同一視しない。

**直し方**: `letI` に**型を明示**して defeq 検査を 1 回だけ走らせる。

```lean
letI : Normal K'.carrier
    (IntermediateField.adjoin K'.carrier ({ᾱ g.gen} : Set K'.closure)) :=
  (stageGeneratorMap α ᾱ hfwd hN g).normal
```

★同じ罠が `FiniteDimensional` でも出る（同じシフトで 2 件）。

## #279 `rw [hc]` は **`hc` の左辺が右辺の中にも現れる**とき自分を書き換える（2026-09-08、同上）

**症状**: `hc : y = c * Φ.symm y * c⁻¹` から `Φ.symm y = c⁻¹ * y * c` を出そうとして
`by rw [hc]; group` と書いた。

```
error: unsolved goals
hc : y = c * Φ.symm y * c⁻¹
⊢ Φ.symm (c * Φ.symm y * c ^ (-1)) = Φ.symm y
```

**原因**: `rw [hc]` は**ゴール全体**の `y` を置き換えるので、左辺の `Φ.symm y` の中の `y` も
置き換わってしまう。★`group` は右辺だけ整理して終わる。

**直し方（★一番安いのは「その補題を使わない」）**: 実測ではこの `have` 自体が要らず、
既に証明済みの一般形（`map_eq_map_of_conj Φ (MulEquiv.refl Γ) H h`）に流し込むだけで済んだ:

```lean
theorem map_eq_self_of_conj (Φ : Γ ≃* Γ) (H : Subgroup Γ) [H.Normal]
    (h : ∀ g : Γ, ∃ c : Γ, Φ g = c * g * c⁻¹) : Subgroup.map (Φ : Γ →* Γ) H = H := by
  rw [map_eq_map_of_conj Φ (MulEquiv.refl Γ) H h]
  ext x
  simp
```

★どうしても要るなら `set z := Φ.symm y with hz` で **左辺を別の fvar にしてから** `rw` する。

## #280 ★★`.cache/mathlib-index.txt` は **`to_additive` の生成名を持たない** —— 加法版の「grep 0 件」は不在の証拠にならない（2026-09-08、pGC 巡回 p 次の跳び）

★**この 1 件で「自前で 60 行書く」が「mathlib を 1 行呼ぶ」に変わった。**

`.cache/mathlib-index.txt` は**ソースに書かれた宣言**から作る。`to_additive` が生成する
加法版は**ソースに名前が出てこない**ので、索引に**1 件も載らない**。

```
grep -n "nnnorm_sum_eq_sup\|norm_sum_eq_sup" .cache/mathlib-index.txt
  → IsUltrametricDist のものは 0 件（BoundedContinuousFunction 等の別物だけ）
grep -n "Analysis/Normed/Group/Ultra.lean" .cache/mathlib-index.txt
  → ★乗法版 IsUltrametricDist.nnnorm_prod_eq_sup_of_pairwise_ne は在る
```

⇒ ★**乗法版が索引に在ったら、加法版の名前を機械的に作って `#check` を投げる**
（`prod → sum`, `mul → add`, `inv → neg`, `div → sub`, `one → zero`）。

```
#check @IsUltrametricDist.nnnorm_sum_eq_sup_of_pairwise_ne
→ ∀ {M ι} [SeminormedAddCommGroup M] [IsUltrametricDist M] {s : Finset ι} {f : ι → M},
    ((↑s).Pairwise fun i j ↦ ‖f i‖₊ ≠ ‖f j‖₊) → ‖∑ i ∈ s, f i‖₊ = s.sup fun i ↦ ‖f i‖₊
```

★**外した名前はエラーで即分かる**ので、`#check` を 3〜4 本まとめて 1 往復で投げてよい:

```
error(lean.unknownIdentifier): Unknown constant `IsUltrametricDist.norm_sum_le_of_forall_le`
```

（正しくは `IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg`。★同じ往復で判明した。）

## #281 `Finset.le_sup` を `exact` で打つと `OrderBot ?m` で止まる —— `apply` にする（2026-09-08、同上）

```
error: typeclass instance problem is stuck
  OrderBot ?m.279
```

`exact Finset.le_sup hjt` は、`sup` の値域（ここでは `ℝ≥0`）を**ゴールと単一化する前に**
インスタンス探索を始めるので詰まる。★`exact` を `apply` に変えるだけで通る。

```
×  · exact Finset.le_sup hjt
○  · apply Finset.le_sup hjt
```

★同じ形は `Finset.sup_le` / `Finset.sup_mono` でも起きうる（`OrderBot` を要求する補題全部）。

## #282 `omit [Inst] in` は **docstring より前**に置く（2026-09-08、同上）

```
error: unexpected token 'omit'; expected 'lemma'
```

★`/-- … -/` と `theorem` の**間**に `omit … in` を挟むと落ちる。docstring は宣言に
直接くっついていなければならない。

```
×  /-- doc -/
   omit [IsUltrametricDist L] in
   theorem foo …

○  omit [IsUltrametricDist L] in
   /-- doc -/
   theorem foo …
```

★そもそもこれが要るのは、`variable` に書いたインスタンスが使われないと
`warning: automatically included section variable(s) unused in theorem …` が出るため。

## #283 `exact_mod_cast` は `‖x‖₊ = ‖y‖₊` → `‖x‖ = ‖y‖` を**渡らない**（2026-09-08、同上）

```
error: mod_cast has type
  ‖x‖₊ = ‖y‖₊
but is expected to have type
  ‖x‖ = ‖y‖
```

★`nnnorm` は「`norm` に coe が付いた形」ではなく**別の関数**なので、`norm_cast` の
単純化では橋が架からない（`coe_nnnorm` は simp 補題であって cast 補題ではない）。

```
×  fun x y h he => h (by exact_mod_cast he)
○  intro x y h he
   exact h (by simpa using congrArg NNReal.toReal he)
```

★**逆向き**（`ℝ` の等式から `ℝ≥0` の等式）は `exact_mod_cast` で通る:

```
have hc : ((‖v j‖₊ * ‖π‖₊ : NNReal) : ℝ) = ((‖u - π‖₊ * ‖w j‖₊ : NNReal) : ℝ) := by
  push_cast
  exact h            -- h は ℝ の等式
exact_mod_cast hc    -- ℝ≥0 の等式が出る
```

★不等式も同様で、`‖x‖₊ ≤ ‖y‖₊` から `‖x‖ ≤ ‖y‖` は
`simpa using NNReal.coe_le_coe.mpr hle` が最短。

## #284 ★★構造体フィールドの型（`A.G` のような射影）に入った値は **関数として適用できない** —— ★**型注釈 `(e : T) x` では直らない**（2026-09-08、pGC Theorem 4.2 の全単射性）

`FilteredGroup.Iso` の `equiv` は `ContinuousMulEquiv A.G B.G` なので `F.equiv g` の型は
`(filtOf RF K').G` になる。これは `K'.absGal`（= `AlgEquiv`）と**定義的に等しい**が、
`x` を渡すとエラボレータが止まる:

```
error: Function expected at
  F.equiv g
but this term has type
  (filtOf RF K').G

Note: Expected a function because this term is being applied to the argument
  x
```

★★**型注釈を足しても、同じエラーが逐語で再現する**。実測（同日）:
`(F.equiv g : K'.absGal) x` と書き換えて再検査したが、
`error: Function expected at / F.equiv g / but this term has type / (filtOf RF K').G` が
**1 文字も変わらずに 2 度目も出た**。★`(e : T)` は `e` を期待型 `T` で
エラボレートするだけで、結果の項の**推論型は据え置かれる**ため、
`CoeFun` の探索は元の型（射影の形）で走る。

```
×  (hF : ∀ g x, F.equiv g x = ᾱ (g (ᾱ.symm x)))
×  (hF : ∀ g x, (F.equiv g : K'.absGal) x = ᾱ (g (ᾱ.symm x)))   -- ★同じエラー
○  (hF : ∀ g, (F.equiv g : K'.absGal) = galMulEquivOf α ᾱ hfwd g)  -- 適用をやめる
○  (show K'.absGal from F.equiv g) x                                -- どうしても点で書くとき
```

★**いちばん安いのは「適用をやめて、群の元の等式で述べる」**（点ごとの形が要る側は
`AlgEquiv.ext` を 1 回打てば渡れる）。★`obtain` で出てきた変数が同じ病気のときは
`set c : K'.absGal := c₀ with hcdef` で**型を貼り替える**と、以後 `c (…)` が通る。

## #285 `AlgEquiv` の coe と `.toRingEquiv` の coe は defeq だが `rw` は跨がない（2026-09-08、同上）

`c : K'.absGal` に対し `c x` と `c.toRingEquiv x` は defeq。しかし `rw` は構文照合なので:

```
error: Tactic `rewrite` failed: Did not find an occurrence of the pattern
  c.toRingEquiv ((extendToClosure α) (τ ((algebraMap K.carrier K.closure) x)))
in the target expression
  c ((extendToClosure α) (τ ((algebraMap K.carrier K.closure) x))) = (algebraMap K'.carrier K'.closure) (β x)
```

★原因は、抽象核を `RingEquiv` で書いて `c.toRingEquiv` を渡したので**結論だけ
`.toRingEquiv` の形**になり、`have` の側が `AlgEquiv` の coe で書かれていたこと。

```
×  have hL : c (… ) = …          -- rw [hτ1 …] が「見つからない」
○  have hL : c.toRingEquiv (…) = …   -- 片側に寄せる
   …
   show c (…) = …                 -- 逆向きに戻したいときは show（defeq なので通る）
   exact (extendToClosure α).injective (c.toRingEquiv.injective (hL.trans hR.symm))
```

★**規則**: 抽象核を `RingEquiv` で書いたなら、具体層でも `.toRingEquiv` に**寄せきる**。
`show` で戻すのは 1 箇所だけにする。

## #286 `generalize hn : deg x = n` の**あとに書く `have` は `deg x` ではなく `n` で書く**（2026-09-08、pGC Ax の塔の降下）

`generalize hn : deg x = n` は**ゴールの中の `deg x` を全部 `n` に置き換える**が、
そのあと自分で書く `have` は `deg x` のままなので**噛み合わなくなる**。実測 2 発:

```
error: linarith failed to find a contradiction
…
this : 1 ≤ ∏ k ∈ Finset.Icc 1 (deg x), c k
a✝ : (∏ k ∈ Finset.Icc 1 n, c k) * ε < ε
⊢ False
failed
```

```
error: Application type mismatch: The argument
  hstep1
has type
  ‖x - x'‖ ≤ (∏ k ∈ Finset.Icc 1 (deg x), c k) * ε
but is expected to have type
  ‖x - x'‖ ≤ (∏ k ∈ Finset.Icc 1 (m + 1), c k) * ε
in the application
  max_le hstep1
```

★`hn : deg x = n` は残っているので `rw [hn] at …` で仮説側を寄せられるが、
★**自分で書く `have` の型は最初から `n`（`obtain ⟨m, rfl⟩ : ∃ m, n = m+1` のあとは `m+1`）で書く**のが速い。

```
×  have hstep1 : ‖x - x'‖ ≤ (∏ k ∈ Finset.Icc 1 (deg x), c k) * ε := …
○  have hstep1 : ‖x - x'‖ ≤ (∏ k ∈ Finset.Icc 1 (m + 1), c k) * ε := …
```

★★**そもそも回避できる**: 予算関数 `F : ℕ → ℝ` を外から渡して結論を `F (deg x) * ε` に
すると、`generalize` 後は `F n * ε` の 1 語で済み、`Finset.Icc` の分解が帰納法の中から消える
（`Found/PGC/AxTowerDecay.lean` の `exists_mem_of_descent_budget` は
`exists_mem_of_descent_prod` より **20 行短い**）。

## #287 `a / b⁻¹ = a * b` の名前は `div_inv_eq` では**ない**（2026-09-08、同上）

```
error(lean.unknownIdentifier): Unknown identifier `div_inv_eq`
```

`ε / (p ^ k)⁻¹` を `p ^ k * ε` にしたいとき:

```
×  rw [zpow_neg, div_inv_eq, mul_comm]
○  rw [zpow_neg, zpow_natCast, div_eq_mul_inv, inv_inv, mul_comm]
```

★`div_eq_mul_inv` → `inv_inv` の 2 段に割るのが確実。`field_simp` でも通るが遅い。
★`zpow_natCast` は `a ^ ((n : ℕ) : ℤ) = a ^ n` で、`Padic.norm_eq_zpow_neg_valuation`
（結論が `ℤ` 冪）から `ℕ` 冪に戻すときに必ず要る。

## #288 `Polynomial.hom_eval₂` を `aeval` の形で `have` に書くと通らない（2026-09-08、pGC 中心化群）

```
error: Type mismatch
  Polynomial.hom_eval₂ q ?m.235 ?m.236 y
has type
  ?m.236 (Polynomial.eval₂ ?m.235 y q) = Polynomial.eval₂ (RingHom.comp ?m.236 ?m.235) (?m.236 y) q
but is expected to have type
  σ ((Polynomial.aeval y) q) = Polynomial.eval₂ (σ.toRingHom.comp (algebraMap K.carrier K.closure)) (σ y) q
```

`σ (q(y)) = (ρ q)(σ y)`（係数を `ρ` でひねる）を出すとき、`hom_eval₂` の結論を
**`aeval` で書いた `have` の型**に当てようとすると上で止まる。`aeval` は
`eval₂ (algebraMap ..)` の**略記だが、統一子はここを潜らない**。

```
×  have h1 : σ (Polynomial.aeval y q) = Polynomial.eval₂ (…) (σ y) q :=
     Polynomial.hom_eval₂ q _ _ y
○  have h1 := Polynomial.hom_eval₂ q (algebraMap K.carrier K.closure) σ.toRingHom y
   rw [← Polynomial.aeval_def, hq, map_zero, hcompose] at h1
```

★型を書かずに受けて **`← Polynomial.aeval_def` で後から `aeval` に畳む**のが確実。
★仕上げは `rw [Polynomial.aeval_def, Polynomial.eval₂_map]`
（`eval₂_map : (p.map f).eval₂ g x = p.eval₂ (g.comp f) x` が「係数をひねった多項式」と
「合成した環準同型」を往復させる唯一の橋）。
★`ρ : A ≃ₐ[k] A` を `Polynomial.map` に渡す coe は **`ρ.toRingEquiv.toRingHom`**。


## #289 `pow_succ` と `pow_succ'` は**向きが逆**（`σ^m * σ` を畳むのは `pow_succ`）

```
error: Tactic `rewrite` failed: Did not find an occurrence of the pattern
  ?a * ?a ^ ?n
in the target expression
  (σ ^ m * σ) • x = σ ^ (m + 1) • x
```

`pow_succ : a ^ (n+1) = a ^ n * a`（右に掛ける）／`pow_succ' : a ^ (n+1) = a * a ^ n`（左に掛ける）。
`Function.iterate_succ_apply : f^[n+1] x = f^[n] (f x)` を使うと出てくる形は **`σ^m * σ`** なので、
畳むのは `← pow_succ` である。`← pow_succ'` は上のエラーで止まる。
★エラー文の「pattern `?a * ?a ^ ?n`」がそのまま `pow_succ'` の右辺なので、
**印字されたパターンを見れば `'` の有無を間違えたと分かる**。

```
×  rw [Function.iterate_succ_apply, ih, ← mul_smul, ← pow_succ']
○  rw [Function.iterate_succ_apply, ih, ← mul_smul, ← pow_succ]
```

## #290 `simpa … using le_refl _` は**項の型が `True` に潰れて**落ちる

```
error: Type mismatch: After simplification, term
  le_refl (‖↑n‖ * ‖σ • x - x‖)
 has type
  True
but is expected to have type
  ‖↑n‖ * ‖σ • x - x‖ ≤ ‖↑n‖ * ‖σ • x - x‖
```

`simpa` は**ゴールだけでなく渡した項の型も** `simp` する。`le_refl _`（や `rfl`、`le_rfl`）を渡すと
その型 `a ≤ a` が `simp` で `True` に潰れ、ゴールと合わなくなる。
`≤` の両辺が等式補題で結ばれているときは `simpa` を使わず **`le_of_eq <等式補題>`** を渡す。

```
×  (by simpa [norm_nsmul_closure] using le_refl (‖(n : K.carrier)‖ * ‖σ • x - x‖))
○  (le_of_eq (norm_nsmul_closure K n _))
```

## #291 `padicNormE.norm_p` は**無い** ——「`‖(p : ℚ_p)‖ = p⁻¹`」は `Padic.norm_p`

```
error(lean.invalidField): Invalid field `norm_p`: The environment does not contain `AbsoluteValue.norm_p`, so it is not possible to project the field `norm_p` from an expression
  padicNormE
of type `AbsoluteValue ℚ_[?m.1] ℚ`
```

`padicNormE` は名前空間ではなく**項**（`AbsoluteValue ℚ_[p] ℚ`）なので、
`padicNormE.norm_p` は `AbsoluteValue.norm_p` へのドット記法に解釈される。
ノルムの補題は名前空間 `Padic` の側にある。

```
×  rw [padicNormE.norm_p]
○  rw [Padic.norm_p]          -- ‖(p : ℚ_[p])‖ = (p : ℝ)⁻¹
```

★「The environment does not contain ... so it is not possible to project the field」が出たら
**「`X.f` の `X` が名前空間ではなく項だった」**の合図である（#117 系）。
★同じ罠: `padicNormE.is_norm` などは在るので「`padicNormE.` は全部だめ」ではない。
`AbsoluteValue` に無い名前だけが落ちる。

## #292 「索引に無い ⇒ mathlib に無い」の反例 6 例目 —— `.cache/mathlib-index.txt` は**下限**

```
grep -n "norm_sum_le_of_forall_le" .cache/mathlib-index.txt   → 0 件
grep -n "sum_range_sub"            .cache/mathlib-index.txt   → 別物 1 件だけ（本命が出ない）
grep -n "sum_Ico_eq_sum_range"     .cache/mathlib-index.txt   → 0 件
```

★`node tools/leanfile.mjs` に `#check` を投げると **3 本とも在る**:

* `IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg : 0 ≤ C → (∀ i ∈ s, ‖f i‖ ≤ C) → ‖∑ i ∈ s, f i‖ ≤ C`
* `Finset.sum_range_sub (f : ℕ → G) (n : ℕ) : ∑ i ∈ range n, (f (i + 1) - f i) = f n - f 0`
* `Finset.sum_Ico_eq_sum_range (f) (m n) : ∑ k ∈ Ico m n, f k = ∑ k ∈ range (n - m), f (m + k)`

☆★**実害**: `Found/PGC/AxLemma.lean` の docstring は
「超距離での Multiset 和の一様上界は mathlib に無い(AxLemma.lean が自前で持つ)」と書き、
`norm_multisetSum_le_of_forall_le` を自前で証明した。
★**Finset 版は在る**（Multiset 版が無いのは本当なので、その判断自体は誤りではない。
だが「無い」と書いた語で引いた次の波は Finset 版も自前で書きかねない）。
★**`.absent` を書く前、あるいは「無い」と docstring に書く前に、必ず `#check` を投げること。**

### ★7 例目 + ★★**綴りが変わっていて 0 件になる形**（2026-09-08、Y / wild 深さの降下）

7 例目は上と同じ `IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg`（索引 0 件・実在）。
★**新しいのはこちら** ——

```
grep -nE "\tSubgroup\.relindex" .cache/mathlib-index.txt   → 1 件（別物）
grep -nE "\tSubgroup\.relIndex" .cache/mathlib-index.txt   → 30 本（relIndex_mul_index など）
```

★**現行 mathlib は `Subgroup.relIndex`（I が大文字）に改名済み**。
小文字で引いて 0 件を見た人は「相対指数の補題は無い」と結論しかねない。
書いてしまうと Lean は

```
error(lean.unknownIdentifier): Unknown constant `Subgroup.relindex`
```

としか言わない（★**改名先を教えてくれない**）。
★**大文字小文字を変えてもう一度引く**、あるいは `grep -i` で当てる
（ただし `grep -i` の飲み込みは #232 を読むこと）。

## #293 `← Nat.cast_one` は **ℕ の引き算の中の `1`** まで書き換えて目標を壊す（2026-09-08、Y / 跳びの上界）

★**#240 と同じ族**（`↑(n-1)` と `↑n - 1` の食い違い）だが、**落ち方が違う**ので別に書く。

`((n - 1 : ℕ) : K) + 1 = (n : K)`（`h1 : n - 1 + 1 = n` は手元にある）を

```
rw [← Nat.cast_one, ← Nat.cast_add, h1]
```

で出そうとすると、★**`← Nat.cast_one` が ℕ 側の `n - 1` の中の `1` にも当たって**
`n - ↑1` になり、次の `← Nat.cast_add` が当たらない:

```
error: Tactic `rewrite` failed: Did not find an occurrence of the pattern
  ↑?m + ↑?n
in the target expression
  ↑(n - ↑1) + 1 = ↑n
```

★目標に出ている `↑(n - ↑1)` が動かぬ証拠である（`1` の上にキャストが付いている）。

**直し方: ℕ の等式を先に作って `simpa` で降ろす。**

```lean
have hx : ((n - 1 + 1 : ℕ) : K) = (n : K) := by rw [h1]
simpa using hx
```

★`push_cast` は `Nat.cast_sub` の副条件（`1 ≤ n`）を出しに行くので、
**`n - 1 + 1 = n` を ℕ のまま `omega` で潰してから 1 度だけキャストする**方が短い。

## #294 `Nat.coprime_succ_self_left` は**無い**（2026-09-08、同上）—— `Nat.coprime_sub_self_left` で作る

```
error(lean.unknownIdentifier): Unknown constant `Nat.coprime_succ_self_left`
```

★`.cache/mathlib-index.txt` を `grep -n "coprime_succ_self\|coprime_pred\|succ_coprime"` しても **0 件**。
在るのは `Nat.coprime_sub_self_left {m n} (h : m ≤ n) : Coprime (n - m) m ↔ Coprime n m` と
その `_right` 版だけ（`Data/Nat/GCD/Basic.lean:178`）。

`Nat.Coprime p (p - 1)`（`0 < p`）はこう作る:

```lean
have h1 : p - (p - 1) = 1 := by omega
rw [← Nat.coprime_sub_self_left (Nat.sub_le p 1), h1]
exact Nat.coprime_one_left _
```

★`m := p - 1` を引く向きで使うのがコツ（`m := 1` で引いても目的の形にならない）。

## #295 `leanfile.mjs` は import 先の **olean** を読む —— さっき自分が足した宣言は `Unknown identifier` になる（2026-09-08、Y / wild 深さの降下）

ファイル A に定理を足し、`leanfile.mjs` で A を検査して `ok` を得たあと、
**A を import する検査用ファイル B** を投げると:

```
ABC3/Found/PGC/_probe.lean:49:4: error(lean.unknownIdentifier): Unknown identifier `exists_natDegree_minpoly_descent_div`
ABC3/Found/PGC/_probe.lean:48:9: error: Tactic `rcases` failed: `x✝ : ?m.424` is not an inductive datatype
```

★**2 行目は 1 行目の後始末で出る偽のエラー**である（`obtain` の右辺が解決できないだけ）。
2 行目だけ見ると「`obtain` の分解が間違っている」と読めるので迷う。

**原因**: `leanfile.mjs` は `lake env lean` を呼ぶだけなので、`import` は **olean** を読む。
A の `.lean` を編集しても olean は古いまま。
★A 自身を投げたときは `.lean` を読むので通る —— **通ったのに import 側で落ちる**のはこれ。

**直し方**: 先に `node tools/build.mjs <A のモジュール名>`（実測 9.6 秒）で olean を作ってから B を投げる。

★同じ理由で、**新しい定理を別ファイルから使う前に 1 回だけ build する**のが最短である
（毎回 build するのではない。A を書き終えた 1 回だけ）。

## #296 体の次数（`wildDepth` など）は **中間体ではなく `MulAction.stabilizer` の指数**で測る —— #59 に触らずに済む（2026-09-08、Y / wild 深さの降下）

`[K(x):K] = deg minpoly_K x` を Galois 群の言葉にするとき、
`K⟮x⟯` を `IntermediateField` として作って `IntermediateField.finrank_eq_fixingSubgroup_index`
を使うと、`K.closure` の中間体 `M` の**さらに中の**中間体になり、
★#59 の「中間体 2 層の `rfl` が kernel を止める」に当たる。

★**`MulAction.stabilizer` なら 1 層も作らない**:

```lean
theorem index_stabilizer_eq_natDegree_minpoly
    {F : Type*} [Field F] {E : Type*} [Field E] [Algebra F E]
    [FiniteDimensional F E] [hN : Normal F E] [Algebra.IsSeparable F E] (y : E) :
    (MulAction.stabilizer (E ≃ₐ[F] E) y).index = (minpoly F y).natDegree
```

★部品は 3 つだけ（**実測 1 往復で通った**）:

* `MulAction.index_stabilizer : (stabilizer G x).index = (orbit G x).ncard`
* `Normal.minpoly_eq_iff_mem_orbit E : minpoly F x = minpoly F y ↔ x ∈ orbit Gal(E/F) y`
  （逆向きは `minpoly.eq_of_irreducible_of_monic` で「根 ⇒ minpoly が一致」）
* `Polynomial.card_rootSet_eq_natDegree (hsep) (hN.splits y)`
  （`ncard → Fintype.card` は `Set.ncard_eq_toFinset_card'` + `Set.toFinset_card`）

★これで「深さが 1 下がる」も群論だけで言える:
`Q ≤ stabilizer y` なら `Subgroup.index_dvd_of_le` で
`deg minpoly y ∣ Q.index`、よって `v_p(deg minpoly y) ≤ v_p(Q.index)`。
★**`IntermediateField` を 1 つも作らずに `AxWildDescent` の 1 段が閉じた**
（`Found/PGC/WildDepthDescent.lean` §5・§6 と `Found/PGC/WildDepthFieldDescent.lean`）。

## #297 `.cache/mathlib-index.txt` の行は **section `variable` を含まない** —— 明示引数を 1 つ落として `Application type mismatch`（2026-09-08、Y / `minpoly` の根 = `σ`-軌道）

```
error: Application type mismatch: The argument
  IsUltrametricDist.norm_natCast_le_one ?m.71
has type
  ∀ (n : ℕ), ‖↑n‖ ≤ 1
but is expected to have type
  ‖↑j‖ ≤ 1
in the application
  le_antisymm (IsUltrametricDist.norm_natCast_le_one ?m.71)
```

```
error: Application type mismatch: The argument
  IsUltrametricDist.norm_intCast_le_one ?m.147
has type
  ∀ (z : ℤ), ‖↑z‖ ≤ 1
but is expected to have type
  ‖↑(j.gcdA p)‖ ≤ 1
```

★**読み方**: 「引数を 1 つ渡したのにまだ `∀` が残っている」＝
★**自分が知らない明示引数が先頭にもう 1 つある**。

★**原因**: `.cache/mathlib-index.txt` は**宣言の字面だけ**を持っていて、
その宣言を囲む `variable (R : Type*) [NormedRing R] [IsUltrametricDist R]` を持たない。
索引の行は

```
lemma IsUltrametricDist.norm_natCast_le_one	Analysis/Normed/Ring/Ultra.lean:67	lemma norm_natCast_le_one (n : ℕ) : ‖(n : R)‖ ≤ 1
```

と出るので `(n : ℕ)` だけに見えるが、★実際の signature は
`norm_natCast_le_one (R) (n : ℕ)` である（`R` が**明示** `variable`）。

★**直し方**: `IsUltrametricDist.norm_natCast_le_one L j` /
`IsUltrametricDist.norm_intCast_le_one L _` のように**担い手の型を先頭に足す**。
★迷ったら全部 `_` にして `exact?` ではなく `#check @名前` で本当の signature を見る
（`@` を付けると section variable も含めて全部出る）。

★★**逆向きの罠も同じ日に踏んだ**: 索引の
`theorem modByMonic_add_div (p q : R[X]) : p %ₘ q + q * (p /ₘ q) = p` は
★**こちらは正しい**（`Monic` の仮定は無く、非 monic では両辺とも `p` で成り立つ）。
`Polynomial.modByMonic_add_div g hm`（`hm : Monic`）と書いて

```
error: Application type mismatch: The argument
  hm
has type
  (minpoly K π).Monic
of sort `Prop` but is expected to have type
  K[X]
of sort `Type u_1`
```

★**索引を疑う前に、まず索引どおりに渡してみる**こと。
「`_ %ₘ _` は monic 仮定が要るはず」という思い込みの方が外れていた。

## #298 自分で `[TopologicalSpace X]` を `variable` に書くと **`krullTopology` を潰す** —— `inst✝` と名前つきインスタンス項の食い違い（2026-09-08、N1 / `ker(Art) ⊓ I_K`）

```
error: Application type mismatch: The argument
  Subgroup.topologicalClosure_minimal H le_rfl hH
has type
  @Subgroup.topologicalClosure Gal(E/k) inst✝¹ AlgEquiv.aut inst✝ H ≤ H
but is expected to have type
  @Subgroup.topologicalClosure Gal(E/k) (krullTopology k E) AlgEquiv.aut ⋯ H ≤ H
```

★**読み方**: 型は同じ `@Subgroup.topologicalClosure Gal(E/k) …` なのに
★**インスタンス引数だけが違う**。片方は `inst✝¹`（＝自分が書いた匿名の
`variable [TopologicalSpace …]`)、もう片方は `(krullTopology k E)`（mathlib の正規のもの）。
★**`inst✝` と「名前のあるインスタンス項」が並んでいたら、ほぼこれである。**

★**原因**: 書いていたのは

```lean
variable {k E : Type*} [Field k] [Field E] [Algebra k E]

theorem foo [IsGalois k E]
    [TopologicalSpace (E ≃ₐ[k] E)] [IsTopologicalGroup (E ≃ₐ[k] E)]   -- ★これが犯人
    {H : Subgroup (E ≃ₐ[k] E)} … := …
```

`E ≃ₐ[k] E` には **`krullTopology k E` という正規のインスタンスが既にある**。
そこへ自前の `[TopologicalSpace _]` を足すと、★**自分の定理の中だけ**別の位相になり、
在庫の補題（`krullTopology` を前提に書かれている）と繋がらなくなる。

★**直し方**: ★**その 2 行を消すだけ**。`IsTopologicalGroup` も
`krullTopology` から自動で見つかる。★1 往復（8〜13 秒）で通った。

★**一般則**: ★**「その型に正規のインスタンスが在るか」を先に測ってから
`variable [C X]` を書く**。`Gal(E/k)`（位相）・`𝒪[K]`（付値位相）・`Units R` は
すべて正規のものが在るので、書いてはいけない。
★逆に `{G : Type*} [Group G] [TopologicalSpace G]` のような**抽象の型変数**なら正しい
（正規のインスタンスが存在しえないので潰しようがない）。

☆★同じ往復でもう 1 つ踏んだのは **#297**（`instIsClosedTopologicalClosureCommutator`
は索引の行に出ない明示引数 `(G)` を持つ）。
★`Subgroup.isClosed_topologicalClosure _` に替えて回避した。

## #299 import 不足の `Nat.Prime.foo` は **`Unknown constant` ではなく `Invalid field … Irreducible.foo`** で出る（2026-09-08、Y / 巡回 `p` 次の層）

★逐語（`hp : p.Prime`、`p : ℕ` に対して `hp.dvd_choose_self` と書いたとき）:

```
Invalid field `dvd_choose_self`: The environment does not contain
`Irreducible.dvd_choose_self`, so it is not possible to project the field
`dvd_choose_self` from an expression
  hp
of type `Irreducible p`
```

★★**「型が `Irreducible p` になっている」と読んではいけない。**
`Nat.Prime p` は定義が `Irreducible p` なので、`Nat.Prime` 名前空間に
その名前が**無いとき**だけ Lean は展開して `Irreducible` 名前空間を探し、
そこにも無いのでこのエラーになる。★**原因は import 不足**（#68 の別の顔）。

★`(hp : Nat.Prime p)` と書き直しても**同じエラーが出る**（型の綴りの問題ではない）。
★直し方: **完全名で書く** —— `Nat.Prime.dvd_choose_self hp` にすると
今度は `Unknown constant `Nat.Prime.dvd_choose_self`` に変わり、import 不足だと分かる。
`grep -n "dvd_choose" .cache/mathlib-index.txt` が
`Data/Nat/Choose/Dvd.lean:35` を指すので `import Mathlib.Data.Nat.Choose.Dvd` で解決した。

★★**同じ形は `Nat.Prime` 以外でも起きる**（`def` で書かれた述語すべて）。
「ドット記法が知らない名前空間を名指してきたら import を疑う」が合図である。

## #300 `to_additive` が生成した名前は `.cache/mathlib-index.txt` に**載らない**（2026-09-08、同上）

★索引を `IsUltrametricDist.norm_` で引くと**乗法版しか出ない**:

```
IsUltrametricDist.norm_pow_le  Analysis/Normed/Group/Ultra.lean:139
  lemma norm_pow_le (x : S) (n : ℕ) : ‖x ^ n‖ ≤ ‖x‖
```

★ところが `#check @IsUltrametricDist.norm_nsmul_le` は**在る**:
`(x : S) (n : ℕ) : ‖n • x‖ ≤ ‖x‖`。
★理由: `decl-index.mjs` は宣言の**字面**しか読まないので、
`@[to_additive nnnorm_nsmul_le]` のように**属性の中にだけ書かれた加法版の名前**を拾えない。

★★**手順**: 加法版が欲しいときは
(i) 索引で**乗法版**（`_pow_` / `_mul_` / `_prod_` / `_div_`）を引き、
(ii) 名前を機械的に読み替えて（`pow→nsmul`, `zpow→zsmul`, `mul→add`, `prod→sum`, `div→sub`）
(iii) **`#check` で確かめる**。
★これで「索引に無い ⇒ 不在」の誤判定を 1 件回避した（`norm_nsmul_le`）。

## #301 `Fact (Nat.Prime 3) := ⟨by norm_num⟩` は **`unsolved goals ⊢ Nat.Prime 3`** で落ちる（2026-09-08、正規化した跡の道）

具体的な素数で `Fact` を立てるとき、`norm_num` が原始性を落とさないことがある。逐語:

```
ABC3/Found/PGC/NormalizedTraceDescent.lean:105:33: error: unsolved goals
⊢ Nat.Prime 3
```

★`norm_num` の素数拡張は `Mathlib.Tactic.NormNum.Prime` に在り、
`ABC3.Found.PGC.AxEpsilonDecay` の import 連鎖には**入っていない**。
★`Unknown constant` ではなく `unsolved goals` で出るので #68 の合図には見えない。

★★**直し方（import を足さずに済む）**: mathlib が持っている名前つきの証明を使う。

```lean
haveI : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩   -- ok
```

`Nat.prime_two` / `Nat.prime_three` / `Nat.prime_five` / `Nat.prime_seven` /
`Nat.prime_eleven` は `Mathlib/Data/Nat/Prime/Defs.lean` に在る（`decide` を書く必要も無い）。
測ったコマンド: `grep -n "Nat.prime_three" .cache/mathlib-index.txt`。

## #302 ℝ の有限積に `Finset.single_le_prod'` は使えない —— `MulLeftMono ℝ`（2026-09-08、多段降下の台帳）

「積の 1 因子は積以下」を `Finset.single_le_prod'` で書くと、ℝ では逐語:

```
ABC3/Found/PGC/WildDescentMultiStep.lean:124:9: error(lean.synthInstanceFailed): failed to synthesize instance of type class
  MulLeftMono ℝ
```

```
ABC3/Found/PGC/WildDescentMultiStep.lean:124:53: error: Application type mismatch: The argument
  zero_le_one
has type
  0 ≤ 1
but is expected to have type
  1 ≤ 1
```

★2 つ目は仮説が `0 ≤ f i` ではなく **`1 ≤ f i`** であることの合図（順序付き**モノイド**の補題だから）。
★1 つ目が本体で、**ℝ は乗法について順序モノイドではない**（負数がある）ので
`Finset.single_le_prod'` / `Finset.prod_le_prod'` 族は ℝ には当たらない。

★★**直し方（`nlinarith` に落とす）** —— 1 因子を `Finset.mul_prod_erase` で外に出す:

```lean
have hmem : d ∈ s := by simp only [Finset.mem_Icc]; omega
have hsplit := Finset.mul_prod_erase s c hmem      -- c d * ∏ (s.erase d) = ∏ s
have h1 : (1:ℝ) ≤ ∏ k ∈ s.erase d, c k := Finset.one_le_prod (fun i _ => hc i)
rw [← hsplit]
nlinarith [hc d]
```

★`Finset.one_le_prod`（`1 ≤ f i` から `1 ≤ ∏`）は ℝ で**そのまま通る**
（順序付き半環の補題で、乗法モノイドの順序を要求しない）。ここが分かれ目。

★同じ命令で名前が 2 つ動いていた（索引の「無いが嘘」ではなく本当に無い）:
`` Unknown constant `Nat.Icc_succ_left` `` → `Finset.Icc 1 d = Finset.Ioc 0 d` は
`by ext k; simp only [Finset.mem_Icc, Finset.mem_Ioc]; omega` で作る
（そのうえで `Finset.prod_Ioc_consecutive` が `∏_{(0,d']}·∏_{(d',d]} = ∏_{(0,d]}` をくれる）。
`` Unknown identifier `le_or_lt` `` → `Nat.lt_or_ge a b : a < b ∨ a ≥ b` を使う。

## #303 `generalize hn : deg x = n` の強帰納では `rw [hn] at *` ではなく `subst hn`（2026-09-08、深い降下の反例）

```
error: Application type mismatch: The argument
  hstep'
has type
  c (deg x) ≤ ∏ k ∈ Finset.Icc (b + 1) (deg x), c k
but is expected to have type
  c (deg x) ≤ ∏ k ∈ Finset.Icc (b + 1) n, c k
in the application
  mul_le_mul_of_nonneg_right hstep'
```

`generalize hn : deg x = n; induction n using Nat.strong_induction_on generalizing x` の枝の中は
ゴールが `n`、`have` で作った補題が `deg x` になり、**両方が同時に見えている**ので
`exact` が上のように落ちる。★`rw [hn] at *` で揃えようとすると `hn` 自身や `ih` まで巻き込んで
別の `Application type mismatch` に化ける（実際に 2 往復とかした）。

★直し方は `intro` を済ませた直後に **`subst hn` の 1 行**。`n` が `deg x` に消え、
`ih : ∀ m < deg x, …` も自動で揃う。以後 `hn` は書かない。

★ついでに: この形の帰納で `hstep` から得た `x₁` は `deg x₁ < deg x` しか言えず、
**`b ≤ deg x₁` は言えない**（一気に閾値の下へ落ちうる）。
`Finset.Icc (b+1) (deg x₁)` を割る補題に `b ≤ deg x₁` を要求すると
`omega could not prove the goal` になるので、`Nat.lt_or_ge (deg x₁) (b+1)` で
**空積の枝を分ける**。空積側は `Finset.Icc_eq_empty` + `omega`。

## #304 `field_simp` は `(1/a)^n` を `(a)⁻¹ ^ n` のまま残す —— 先に `one_div_pow`（2026-09-08、対の予算の閉じた形）

`axDecay p i = p ^ ((1/(p-1)) * (1/p)^(i-1))` の有限積を閉じた形にする帰納で、
`congr 1` のあと `field_simp; ring` を撃つと **`ring` が閉じない**:

```
error: unsolved goals
case succ.e_a
⊢ -(↑p * ↑p ^ (d * 2) * ↑p ^ (k' * 2) * (↑p)⁻¹ ^ d * (↑p)⁻¹ ^ k') - ↑p ^ 2 * ↑p ^ d * ↑p ^ k' +
        ↑p ^ 2 * ↑p ^ (d * 2) * ↑p ^ k' +
Try this:
  [apply] ring_nf
  
  The `ring` tactic failed to close the goal. Use `ring_nf` to obtain a normal form.
```

★読み方: 正規形に **`↑p ^ (d*2)`（分母を払った跡）と `(↑p)⁻¹ ^ d`（払えなかった跡）が同居**している。
`field_simp` は `1/a` は分母として拾うが、**`(1/a)^n`（`n` が変数）は拾わない**。

★直し方は `field_simp` の**前**に `one_div_pow` を 1 個足すだけ:

```lean
rw [..., ← Real.rpow_add (by positivity), one_div_pow]
congr 1
field_simp
ring
```

`one_div_pow : (1 / a) ^ n = 1 / a ^ n`。これで分母が `a ^ n` の 1 個になり `field_simp` が払える。
★`(1/a)^n` が指数 `n` を**リテラルでなく変数**で持つときだけ起きる（リテラルなら `norm_num` が潰す）。

## #305 `deg` で強帰納する降下核は、結論に **`deg x' ≤ deg x` を足しておく**（2026-09-08、対の予算）

「1 段で `deg` が `d → d'` に落ち、損失は `∏_{k ∈ [d'+1, d]} c k`」という核で、
再帰の戻り値 `x'` が `deg x₁` より**上**に居る場合を潰せず、予算の分割
`(∏_{[deg x'+1, deg x₁]}) * (∏_{[deg x₁+1, deg x]}) = ∏_{[deg x'+1, deg x]}`
の仮定 `deg x' ≤ deg x₁` が出せない。無理に出そうとすると:

```
error: Application type mismatch: The argument
  hcon
has type
  deg x₁ < deg x'
but is expected to have type
  deg x' < deg x₁
in the application
  lt_of_le_of_lt (le_refl (deg x')) hcon
```

★直し方: **結論を 1 つ強くする**。`∃ x', deg x' ≤ b ∧ deg x' ≤ deg x ∧ ‖x - x'‖ ≤ …`
と書けば、基底枝は `le_refl _`、帰納枝は `le_trans hle (le_of_lt h1)` で通り、
再帰の戻り値に `deg x' ≤ deg x₁` が付いてくる。

★★もう 1 つの罠（こちらは**主張が偽になる**）: 着地の予算を
`∏_{k ∈ [b+1, deg x]} c k`（閾値 `b` で書く）とすると**偽**である。
`b = 1`, `deg x = 2` で途中 `deg x' = 0` まで落ちると、支払いは `∏_{[1,2]}` で
許容 `∏_{[2,2]}` を超える。★正しくは**着地の `deg x'` で書く**:
`∏_{k ∈ [deg x' + 1, deg x]} c k`。

## #306 scratch の Lean に `import Mathlib` と書くと `leanfile.mjs` が 354 秒（プロジェクトの module なら 8.6 秒）

試作を scratchpad の `.lean` に置いて `node tools/leanfile.mjs <file>` で回すとき、
先頭を `import Mathlib` にすると 1 往復が **354.5 秒**になる（2026-09-08 実測）。

```
ok  .../scratchpad/grad/T4.lean  —— 354.5 秒
```

そのうえ

```
Command did not complete within its 120s timeout and was moved to the background
```

で背景化するので、**待ちの往復がさらに 2 回増える**（同じ試作を
`import ABC3.Found.PGC.DeepDescentPairDirect` に変えたら **8.9 秒**、40 倍差）。

★直し方: 試作でも**その持ち場が最終的に import するプロジェクトの module** を import する。
ABC3 の module は Mathlib を推移的に持つので、`Mathlib` 全体の olean を読む理由はまず無い。
★`--similar` で先行節を探すときは **`--similar <下書きのファイル>`**（文字列ではなくパス）:
文字列を渡すと `Error: ENOENT: no such file or directory, open 'D:\Math_ABC3\unsolved goals'`。

## #307 `positivity` は `↑p - 1` の引き算で止まる／ℕ の台帳から正の因子を約す 3 行（2026-09-08、跳びと different のトレードオフ）

**失敗形**（`Found/PGC/JumpDefectTradeoff.lean`、`2 ≤ p` は**仮定**にある）

```lean
have hpos : (0:ℤ) < ((p:ℤ) - 1) * (p:ℤ) ^ 2 := by positivity
```

```
error: failed to prove positivity/nonnegativity/nonzeroness
```

★`positivity` は**式の形だけ**を見るので、仮定 `2 ≤ p` から来る `0 < ↑p - 1` は見えない
（先行節「`positivity` が `≠ 0` を出せないとき」と同じ根。あちらはノルム、こちらは引き算）。

**直し方**: 引き算の因子だけ `linarith` に回し、残りは `positivity` に任せる。

```lean
have hpos : (0:ℤ) < ((p:ℤ) - 1) * (p:ℤ) ^ 2 := mul_pos (by linarith) (by positivity)
```

★ついでに（同じ証明で必要になった)。ℕ の台帳
`(p-1)^2 * p^(0+2) * J ≤ (p^2-1) * p * (p^2*e)` から正の因子 `(p−1)p²` を約すのは

```lean
zify [h1, h2]                       -- h1 : 1 ≤ p, h2 : 1 ≤ p ^ 2（ℕ の引き算を外す）
have e1 : LHS = (X) * (((p:ℤ)-1) * (p:ℤ)^2) := by ring   -- ★右に因子を寄せる
have e2 : RHS = (Y) * (((p:ℤ)-1) * (p:ℤ)^2) := by ring
rw [e1, e2, mul_le_mul_iff_left₀ hpos]
```

の 3 行。★`nlinarith` に丸投げすると `linarith failed to find a contradiction` で落ちる
（約分は非線形なので、`↑p ^ (0 + 2)` が残っていると特にだめ。`ring` を通す `show`/`have` で
指数の `0 + 2` ごと潰すのが早い）。★名前は `mul_le_mul_left` ではない——それは
`(bc : b ≤ c) (a : α) : b * a ≤ c * a` で iff ではない。`mul_le_mul_iff_left₀ (a0 : 0 < a) : b * a ≤ c * a ↔ b ≤ c`
を `grep -n "mul_le_mul_iff" .cache/mathlib-index.txt`（0.2 秒）で引く。
