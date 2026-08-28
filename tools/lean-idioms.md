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
